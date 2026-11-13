#!/usr/bin/env bash
# deplete_host.sh — stream an aligned BAM/CRAM to FASTQ and remove host reads with deacon.
#
# All reads are extracted (not just unmapped ones) so depletion does not depend on the
# original alignment. Mates are re-paired with `samtools collate`, written as an
# interleaved stream into `deacon filter --deplete`, and split into R1/R2 outputs.
# Reads whose mate is absent (singletons) and unpaired reads are depleted separately.
#
# Outputs (in OUTDIR):
#   <prefix>.host_depleted.R1.fq.gz
#   <prefix>.host_depleted.R2.fq.gz
#   <prefix>.host_depleted.unpaired.fq.gz
#   <prefix>.deacon.paired.json      deacon summary (paired)
#   <prefix>.deacon.unpaired.json    deacon summary (singletons + unpaired)
#   <prefix>.depletion_stats.tsv     one-line summary for aggregation
#
# Any failure anywhere in the pipe (e.g. a truncated/corrupt BAM) makes this script exit non-zero.
set -euo pipefail

usage() {
  cat >&2 <<EOF
Usage: $(basename "$0") -i INPUT.(bam|cram) -x INDEX.idx -p PREFIX [options]
  -i  input BAM or CRAM (coordinate- or name-sorted)
  -x  deacon index (e.g. panhuman-1.k31w15.idx)
  -p  output prefix (sample name)
  -o  output directory                      [.]
  -r  reference FASTA (required for CRAM)
  -t  threads                               [nproc]
  -a  deacon --abs-threshold                [2]
  -R  deacon --rel-threshold                [0.01]
  -T  temp directory for samtools collate   [OUTDIR/tmp]
EOF
  exit 2
}

INPUT="" INDEX="" PREFIX="" OUTDIR="." REF="" THREADS="$(nproc)" ABS=2 REL=0.01 TMPD=""
while getopts "i:x:p:o:r:t:a:R:T:h" opt; do
  case $opt in
    i) INPUT=$OPTARG ;; x) INDEX=$OPTARG ;; p) PREFIX=$OPTARG ;; o) OUTDIR=$OPTARG ;;
    r) REF=$OPTARG ;; t) THREADS=$OPTARG ;; a) ABS=$OPTARG ;; R) REL=$OPTARG ;;
    T) TMPD=$OPTARG ;; *) usage ;;
  esac
done
[[ -n $INPUT && -n $INDEX && -n $PREFIX ]] || usage
[[ -f $INPUT ]] || { echo "ERROR: input not found: $INPUT" >&2; exit 1; }
[[ -f $INDEX ]] || { echo "ERROR: index not found: $INDEX" >&2; exit 1; }
if [[ $INPUT == *.cram && -z $REF ]]; then
  echo "ERROR: CRAM input requires -r reference FASTA" >&2; exit 1
fi

mkdir -p "$OUTDIR"
TMPD=${TMPD:-$OUTDIR/tmp}
mkdir -p "$TMPD"
REF_OPT=()
[[ -n $REF ]] && REF_OPT=(--reference "$REF")

out="$OUTDIR/$PREFIX.host_depleted"
log() { echo "[$(date '+%F %T')] $*" >&2; }

for t in samtools deacon; do command -v $t >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }; done
log "samtools $(samtools --version | head -1 | cut -d' ' -f2), $(deacon --version)"
log "input=$INPUT index=$INDEX threads=$THREADS abs=$ABS rel=$REL"

# Fail fast on obviously broken inputs (missing EOF block etc.); -u allows unaligned BAMs (no @SQ)
samtools quickcheck -u -v "$INPUT" || { echo "ERROR: samtools quickcheck failed for $INPUT" >&2; exit 1; }

# 1) collate mates (fast mode: primary alignments only), 2) FASTQ, 3) deacon deplete.
#    samtools fastq: paired READ1/READ2 -> stdout (interleaved, /1 /2 suffixes),
#    singletons -> -s, reads flagged neither/both READ1/READ2 -> -0.
#    Default -F 0x900 drops secondary/supplementary so each read appears once; duplicates
#    and QC-fail reads are kept on purpose (microbial reads are often flagged as dups).
log "extracting + depleting paired reads"
samtools collate -@ "$THREADS" -f -r 100000 -O -u "${REF_OPT[@]}" "$INPUT" "$TMPD/collate" \
  | samtools fastq -@ "$THREADS" "${REF_OPT[@]}" \
      -s "$TMPD/singletons.fq.gz" -0 "$TMPD/other.fq.gz" - \
  | deacon filter --deplete --interleaved --check-pairs \
      -a "$ABS" -r "$REL" -t "$THREADS" \
      -s "$OUTDIR/$PREFIX.deacon.paired.json" \
      "$INDEX" - -o "$out.R1.fq.gz" -O "$out.R2.fq.gz"

log "depleting singleton/unpaired reads"
# deacon does not decompress stdin, so decompress here
zcat "$TMPD/singletons.fq.gz" "$TMPD/other.fq.gz" \
  | deacon filter --deplete -a "$ABS" -r "$REL" -t "$THREADS" \
      -s "$OUTDIR/$PREFIX.deacon.unpaired.json" \
      "$INDEX" - -o "$out.unpaired.fq.gz"

rm -rf "$TMPD"

# One-line stats table (reads, not pairs; paired counts include both mates)
python3 - "$PREFIX" "$OUTDIR/$PREFIX.deacon.paired.json" "$OUTDIR/$PREFIX.deacon.unpaired.json" \
  > "$OUTDIR/$PREFIX.depletion_stats.tsv" <<'PY'
import json, sys
name, pj, uj = sys.argv[1:]
p, u = json.load(open(pj)), json.load(open(uj))
def g(d, *keys):
    for k in keys:
        if k in d: return d[k]
    return "NA"
cols = ["sample", "paired_in", "paired_out", "unpaired_in", "unpaired_out", "frac_kept"]
pi, po = g(p, "seqs_in"), g(p, "seqs_out")
ui, uo = g(u, "seqs_in"), g(u, "seqs_out")
tot_in = sum(x for x in (pi, ui) if isinstance(x, int))
tot_out = sum(x for x in (po, uo) if isinstance(x, int))
print("\t".join(cols))
print("\t".join(map(str, [name, pi, po, ui, uo, f"{tot_out / tot_in:.6g}" if tot_in else "NA"])))
PY
cat "$OUTDIR/$PREFIX.depletion_stats.tsv" >&2
log "done"
