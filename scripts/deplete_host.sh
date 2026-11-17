#!/usr/bin/env bash
# deplete_host.sh — stream an aligned BAM/CRAM to FASTQ and remove host reads with deacon.
#
# A read pair is removed if EITHER
#   (1) the original alignment placed both mates on the host reference (alignment filter), OR
#   (2) deacon matches it to the host minimizer index (panhuman-1: HPRC + CHM13 + GRCh38).
# (1) catches repetitive / low-complexity human reads that carry too few distinctive minimizers
# for deacon; (2) catches human sequence absent from the alignment reference. A pair is kept by
# (1) if at least one mate was unmapped; both mates share that decision, so pairs stay intact.
# Disable (1) with -A to deplete on deacon alone.
#
# Surviving primary reads are re-paired with `samtools collate`, streamed interleaved into
# `deacon filter --deplete`, and split into R1/R2. Singletons and unpaired reads are depleted
# separately.
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
  -A  disable the alignment filter (deplete with deacon only)
  -N  comma-separated NON-host contigs: reads mapped to these (or whose mate is) count as
      unmapped for the alignment filter      [$NONHOST]
EOF
  exit 2
}

# Non-human sequences bundled in common human references (GRCh38 analysis sets carry chrEBV).
NONHOST="chrEBV,EBV,NC_007605,NC_007605.1,chrphiX,phiX,phiX174,NC_001422,NC_001422.1"
INPUT="" INDEX="" PREFIX="" OUTDIR="." REF="" THREADS="$(nproc)" ABS=2 REL=0.01 TMPD="" ALN_FILTER=1
while getopts "i:x:p:o:r:t:a:R:T:AN:h" opt; do
  case $opt in
    i) INPUT=$OPTARG ;; x) INDEX=$OPTARG ;; p) PREFIX=$OPTARG ;; o) OUTDIR=$OPTARG ;;
    r) REF=$OPTARG ;; t) THREADS=$OPTARG ;; a) ABS=$OPTARG ;; R) REL=$OPTARG ;;
    T) TMPD=$OPTARG ;; A) ALN_FILTER=0 ;; N) NONHOST=$OPTARG ;; *) usage ;;
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
log "input=$INPUT index=$INDEX threads=$THREADS abs=$ABS rel=$REL alignment_filter=$ALN_FILTER"

# Fail fast on obviously broken inputs (missing EOF block etc.); -u allows unaligned BAMs (no @SQ)
samtools quickcheck -u -v "$INPUT" || { echo "ERROR: samtools quickcheck failed for $INPUT" >&2; exit 1; }

# 0) primary records only (-F 0x900: no secondary/supplementary, so each read appears once);
#    duplicates and QC-fail reads are kept on purpose (microbial reads are often flagged as dups).
#    Counts of primary reads and of reads dropped by the alignment filter go to JSON.
# 1) alignment filter: keep a record only if it or its mate is unmapped.
# 2) collate mates, 3) FASTQ, 4) deacon deplete.
#    Not `collate -f`: fast mode silently drops records without the paired flag (single-end reads
#    that some CBTN CRAMs carry), which the read-accounting check below caught (up to 1.2% of the
#    reads passing the alignment filter). Collate only sees filtered reads, so default mode is cheap.
#    samtools fastq: paired READ1/READ2 -> stdout (interleaved, /1 /2 suffixes),
#    singletons -> -s, reads flagged neither/both READ1/READ2 -> -0.
#    Reads (or mates) on non-host contigs (e.g. chrEBV) are treated as unmapped, so viral reads
#    the aligner placed on them survive to deacon. Both mates evaluate identically
#    (rname/rnext swap), so pairs stay intact.
#    Only non-host contigs present in the header go into the expression: each string comparison
#    is evaluated for every mapped record on a single thread (18 comparisons cut throughput from
#    ~5M to ~0.3M records/s, which bottlenecked whole-genome runs).
ALN_EXPR='1'
if [[ $ALN_FILTER == 1 ]]; then
  ALN_EXPR='flag.unmap || flag.munmap'
  hdr_contigs=$(samtools view -H "$INPUT" | awk -F'\t' '$1=="@SQ"{for(i=2;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}')
  IFS=',' read -ra _nh <<< "$NONHOST"
  for c in "${_nh[@]}"; do
    [[ $c =~ ^[A-Za-z0-9._|-]+$ ]] || { echo "ERROR: bad contig name in -N: $c" >&2; exit 1; }
    grep -qxF -- "$c" <<< "$hdr_contigs" || continue
    ALN_EXPR+=" || rname == \"$c\" || rnext == \"$c\""
  done
fi
log "alignment filter keeps: $ALN_EXPR"
log "extracting + depleting paired reads"
samtools view -@ "$THREADS" -u -F 0x900 "${REF_OPT[@]}" \
    --save-counts "$OUTDIR/$PREFIX.counts.primary.json" "$INPUT" \
  | samtools view -u -e "$ALN_EXPR" --save-counts "$OUTDIR/$PREFIX.counts.alignment_filter.json" - \
  | samtools collate -@ "$THREADS" -O -u - "$TMPD/collate" \
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

# One-line stats table, in reads (paired counts include both mates):
#   primary_reads        primary records in the input
#   aln_removed          removed by the alignment filter (both mates mapped)
#   paired_in/out        into / out of deacon (paired);  unpaired_in/out likewise
#   frac_kept            (paired_out + unpaired_out) / primary_reads
python3 - "$PREFIX" "$OUTDIR/$PREFIX" > "$OUTDIR/$PREFIX.depletion_stats.tsv" <<'PY'
import json, sys
name, base = sys.argv[1:]
ld = lambda s: json.load(open(f"{base}.{s}.json"))
prim, aln = ld("counts.primary"), ld("counts.alignment_filter")
p, u = ld("deacon.paired"), ld("deacon.unpaired")
primary = prim["records_filter_accepted"]
aln_removed = aln["records_filter_rejected"]
pi, po, ui, uo = p["seqs_in"], p["seqs_out"], u["seqs_in"], u["seqs_out"]
# every primary read must be accounted for: alignment-filtered, or sent to deacon
assert aln["records_processed"] == primary, "alignment-filter input != primary reads"
assert aln_removed + pi + ui == primary, (
    f"read accounting mismatch: {aln_removed}+{pi}+{ui} != {primary}")
cols = ["sample", "primary_reads", "aln_removed", "paired_in", "paired_out",
        "unpaired_in", "unpaired_out", "frac_kept"]
row = [name, primary, aln_removed, pi, po, ui, uo,
       f"{(po + uo) / primary:.6g}" if primary else "NA"]
print("\t".join(cols)); print("\t".join(map(str, row)))
PY
cat "$OUTDIR/$PREFIX.depletion_stats.tsv" >&2
log "done"
