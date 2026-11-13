cwlVersion: v1.2
class: CommandLineTool
label: deacon host depletion (BAM/CRAM -> host-depleted FASTQ)
doc: |
  Streams every read of an aligned BAM/CRAM to FASTQ (samtools collate + fastq), re-pairs mates,
  and removes host reads with `deacon filter --deplete` against a minimizer index
  (panhuman-1.k31w15.idx recommended). Outputs paired R1/R2, an unpaired file (singletons),
  deacon JSON summaries and a one-line stats TSV.

  Deploy: sbpack cavatica jkirkland/<project>/deacon-host-depletion cavatica/deacon-host-depletion.cwl

$namespaces:
  sbg: https://sevenbridges.com

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: ghcr.io/jbk708/deacon-host-depletion:0.2.0
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 32000

hints:
  - class: sbg:AWSInstanceType
    # 16 vCPU / 32 GB, 1 TB EBS: room for a ~100 GB CRAM plus collate temp files
    value: c5.4xlarge;ebs-gp2;1024
  - class: sbg:maxNumberOfParallelInstances
    value: 1

inputs:
  input_reads:
    type: File
    doc: Aligned BAM or CRAM (index not required; reads are streamed in full).
    secondaryFiles:
      - pattern: .crai
        required: false
      - pattern: .bai
        required: false
  reference:
    type: File?
    doc: Reference FASTA the CRAM was encoded against (Kids First GRCh38 Homo_sapiens_assembly38.fasta). Required for CRAM.
    secondaryFiles:
      - pattern: .fai
        required: false
  deacon_index:
    type: File
    doc: deacon minimizer index, e.g. panhuman-1.k31w15.idx (index format v3; deacon >= 0.11).
  output_prefix:
    type: string?
    doc: Output prefix; defaults to the input file's basename without extension.
  threads:
    type: int
    default: 16
  abs_threshold:
    type: int
    default: 2
  rel_threshold:
    type: float
    default: 0.01

arguments:
  - shellQuote: false
    valueFrom: >-
      deplete_host.sh
      -i $(inputs.input_reads.path)
      -x $(inputs.deacon_index.path)
      -p $(inputs.output_prefix ? inputs.output_prefix : inputs.input_reads.nameroot)
      -o .
      -t $(inputs.threads)
      -a $(inputs.abs_threshold)
      -R $(inputs.rel_threshold)
      $(inputs.reference ? "-r " + inputs.reference.path : "")
      -T ./tmp

outputs:
  r1:
    type: File
    outputBinding: {glob: "*.host_depleted.R1.fq.gz"}
  r2:
    type: File
    outputBinding: {glob: "*.host_depleted.R2.fq.gz"}
  unpaired:
    type: File
    outputBinding: {glob: "*.host_depleted.unpaired.fq.gz"}
  deacon_summaries:
    type: File[]
    outputBinding: {glob: "*.json"}
  stats:
    type: File
    outputBinding: {glob: "*.depletion_stats.tsv"}
