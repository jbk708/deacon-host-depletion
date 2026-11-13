#!/usr/bin/env bash
# Self-contained end-to-end test with synthetic data (no patient reads).
# Builds a random "host" and "microbe" genome, simulates paired reads from both, packs them into a
# coordinate-style BAM (pairs + singletons), indexes the host with deacon, runs deplete_host.sh,
# and asserts host reads are removed, microbe reads are kept, and R1/R2 stay paired.
# Also checks that a truncated BAM makes the pipeline fail.
# Requires samtools, deacon, python3 and deplete_host.sh on PATH.
set -euo pipefail
work=$(mktemp -d); trap 'rm -rf "$work"' EXIT; cd "$work"

python3 - <<'PY'
import random
random.seed(1)
comp = str.maketrans("ACGT", "TGCA")
def genome(n): return "".join(random.choice("ACGT") for _ in range(n))
host, microbe = genome(300_000), genome(60_000)
open("host.fa", "w").write(">host\n" + host + "\n")
def pairs(g, tag, n, L=150, ins=350):
    out = []
    for i in range(n):
        s = random.randint(0, len(g) - ins)
        f = g[s:s + ins]
        out.append((f"{tag}_{i}", f[:L], f[-L:].translate(comp)[::-1]))
    return out
reads = pairs(host, "HOST", 5000) + pairs(microbe, "MICROBE", 2000)
q = "I" * 150
with open("r.sam", "w") as fh:
    fh.write("@HD\tVN:1.6\tSO:unsorted\n")
    for name, r1, r2 in reads:
        fh.write(f"{name}\t77\t*\t0\t0\t*\t*\t0\t0\t{r1}\t{q}\n")
        fh.write(f"{name}\t141\t*\t0\t0\t*\t*\t0\t0\t{r2}\t{q}\n")
    # singletons: R1 present, mate absent
    for name, r1, _ in pairs(host, "HOSTSINGLE", 300) + pairs(microbe, "MICROBESINGLE", 100):
        fh.write(f"{name}\t73\t*\t0\t0\t*\t*\t0\t0\t{r1}\t{q}\n")
PY
samtools view -b -o test.bam r.sam
deacon index build host.fa > host.idx 2>/dev/null

deplete_host.sh -i test.bam -x host.idx -p t -o out -t 4 2> run.log

r1=out/t.host_depleted.R1.fq.gz r2=out/t.host_depleted.R2.fq.gz un=out/t.host_depleted.unpaired.fq.gz
count() { zcat "$1" | awk -v p="$2" 'NR%4==1 && $1 ~ "^@"p"_" {n++} END{print n+0}'; }
fail=0
check() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fail=1; fi; }
h=$(count $r1 HOST); m=$(count $r1 MICROBE); hs=$(count $un HOSTSINGLE); ms=$(count $un MICROBESINGLE)
echo "kept pairs: host=$h/5000 microbe=$m/2000 | kept singletons: host=$hs/300 microbe=$ms/100"
check "host pairs removed (<=0.5% kept)"      "[ $h -le 25 ]"
check "microbe pairs kept (>=99%)"            "[ $m -ge 1980 ]"
check "host singletons removed"               "[ $hs -le 2 ]"
check "microbe singletons kept"               "[ $ms -ge 99 ]"
mism=$(paste <(zcat $r1 | awk 'NR%4==1{sub(/\/1$/,"",$1);print $1}') <(zcat $r2 | awk 'NR%4==1{sub(/\/2$/,"",$1);print $1}') | awk '$1!=$2' | wc -l)
check "R1/R2 names paired"                    "[ $mism -eq 0 ]"
check "stats TSV written"                     "[ \$(wc -l < out/t.depletion_stats.tsv) -eq 2 ]"

head -c $(( $(stat -c %s test.bam) / 2 )) test.bam > trunc.bam
if deplete_host.sh -i trunc.bam -x host.idx -p tr -o out_tr -t 2 2>/dev/null; then
  echo "FAIL truncated BAM should fail"; fail=1
else
  echo "PASS truncated BAM rejected"
fi
exit $fail
