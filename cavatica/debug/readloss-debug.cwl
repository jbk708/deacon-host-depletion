cwlVersion: v1.2
class: CommandLineTool
label: deacon host depletion - read-loss diagnostic
doc: |
  One-off diagnostic for "read accounting mismatch" failures. Re-runs the pipeline's alignment filter,
  samtools collate and samtools fastq on one CRAM/BAM and reports, per stage, how many reads come out,
  which read names go missing, and their SAM records. No deacon step; outputs are small text files.

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: ghcr.io/jbk708/deacon-host-depletion:0.2.2
  ResourceRequirement:
    coresMin: 16
    ramMin: 32000
  InitialWorkDirRequirement:
    listing:
      - entryname: readloss_debug.sh
        entry: |
          set -euo pipefail
          IN=$1; REF=$2; P=$3; T=16
          EXPR='flag.unmap || flag.munmap || rname == "chrEBV" || rnext == "chrEBV"'
          mkdir -p tmp
          log() { echo "[\$(date '+%F %T')] $*" | tee -a $P.debug_log.txt >&2; }
          log "samtools \$(samtools --version | head -1)"
          # stage 1: primary records -> alignment filter (compressed, small)
          samtools view -@ $T -u -F 0x900 --reference "$REF" --save-counts $P.debug_counts.primary.json "$IN" \
            | samtools view -@ 4 -b -e "$EXPR" --save-counts $P.debug_counts.filter.json -o filt.bam -
          log "filter output records: \$(samtools view -c filt.bam)"
          # stage 2: collate exactly as the pipeline does
          samtools collate -@ $T -f -r 100000 -O -u filt.bam tmp/c > coll.bam
          log "collate output records: \$(samtools view -c coll.bam)"
          # stage 3: fastq exactly as the pipeline does
          samtools fastq -@ $T -s single.fq -0 other.fq coll.bam > paired.fq 2> $P.debug_fastq_stderr.txt
          log "fastq records: paired=\$(( \$(wc -l < paired.fq) / 4 )) single=\$(( \$(wc -l < single.fq) / 4 )) other=\$(( \$(wc -l < other.fq) / 4 ))"
          # per-name record counts at each stage; names whose counts differ are the lost reads
          samtools view filt.bam | cut -f1 | LC_ALL=C sort -S 4G -T tmp | uniq -c | awk '{print $2"\t"$1}' > filt.cnt
          samtools view coll.bam | cut -f1 | LC_ALL=C sort -S 4G -T tmp | uniq -c | awk '{print $2"\t"$1}' > coll.cnt
          cat paired.fq single.fq other.fq | awk 'NR%4==1{n=substr($1,2); sub(/\/[12]$/,"",n); print n}' \
            | LC_ALL=C sort -S 4G -T tmp | uniq -c | awk '{print $2"\t"$1}' > fq.cnt
          LC_ALL=C join -t $'\t' -a1 -e 0 -o 0,1.2,2.2 filt.cnt coll.cnt | awk -F'\t' '$2!=$3' > $P.debug_lost_in_collate.tsv
          LC_ALL=C join -t $'\t' -a1 -e 0 -o 0,1.2,2.2 coll.cnt fq.cnt | awk -F'\t' '$2!=$3' > $P.debug_lost_in_fastq.tsv
          log "names with lost records: collate=\$(wc -l < $P.debug_lost_in_collate.tsv) fastq=\$(wc -l < $P.debug_lost_in_fastq.tsv)"
          # records (from the filter output) of every affected name, plus a flag histogram
          cut -f1 $P.debug_lost_in_collate.tsv $P.debug_lost_in_fastq.tsv | sort -u > lost.names
          if [ -s lost.names ]; then
            samtools view -N lost.names filt.bam > $P.debug_lost_records.sam
            cut -f2 $P.debug_lost_records.sam | sort | uniq -c | sort -rn > $P.debug_lost_flag_hist.txt
            # how many records per affected name in the filter output (pairs vs. >2 records per name)
            cut -f1 $P.debug_lost_records.sam | sort | uniq -c | awk '{print $1}' | sort -n | uniq -c > $P.debug_lost_records_per_name.txt
          fi
          head -c 2000000 $P.debug_lost_records.sam > $P.debug_lost_records.head.sam 2>/dev/null || true
          rm -f $P.debug_lost_records.sam
          log done

inputs:
  input_reads:
    type: File
    secondaryFiles:
      - pattern: .crai
        required: false
  reference:
    type: File
    secondaryFiles:
      - pattern: .fai
        required: false
  output_prefix:
    type: string

arguments:
  - shellQuote: false
    valueFrom: >-
      bash readloss_debug.sh $(inputs.input_reads.path) $(inputs.reference.path) $(inputs.output_prefix)

hints:
  - class: sbg:AWSInstanceType
    value: c5.4xlarge;ebs-gp2;256

$namespaces:
  sbg: https://sevenbridges.com

outputs:
  debug_files:
    type: File[]
    outputBinding: {glob: "*.debug_*"}
