# deacon-host-depletion

Host (human) read depletion for aligned WGS BAM/CRAM files using
[deacon](https://github.com/bede/deacon), packaged for two platforms:

| Platform | Data | Wrapper |
|---|---|---|
| DNAnexus (St. Jude Cloud) | St. Jude WGS BAMs (`azure:westus`) | `dnanexus/deacon_host_depletion/` applet |
| Cavatica (Seven Bridges) | CBTN / PBTA WGS CRAMs | `cavatica/deacon-host-depletion.cwl` + Docker image |

Both call the same script, `scripts/deplete_host.sh`.

## What it does

A read pair is removed if **either** test says it is human:

1. **Alignment filter:** the original alignment placed both mates on the reference
   (`samtools view -e 'flag.unmap || flag.munmap'` keeps pairs with at least one unmapped mate).
   This catches repetitive and low-complexity human reads that carry too few distinctive
   minimizers for deacon.
2. **deacon:** the pair matches the panhuman-1 minimizer index. This catches human sequence
   missing from the alignment reference (pangenome / T2T regions).

Steps: primary records (duplicates and QC-fail reads kept) → alignment filter →
`samtools collate -f | samtools fastq` → `deacon filter --deplete --interleaved --check-pairs` →
R1/R2 with mates kept in sync. Singletons (mate missing) and unpaired reads go through deacon
separately. Use `-A` to skip the alignment filter and deplete on deacon alone.

**Non-host contigs.** Some human references carry non-human sequence (the Kids First
`Homo_sapiens_assembly38.fasta` used for CBTN CRAMs includes `chrEBV`). Reads, or mates of reads,
mapped to these are treated as unmapped by the alignment filter, so they reach deacon instead of
being discarded. Default: EBV and phiX contig names; override with `-N chrEBV,...`. The St. Jude
`hg38m1x` reference has no non-human contigs.

On a pilot germline WGS sample (SJAML030881_G1), deacon alone kept 1.05M read pairs the original
alignment had placed on hg38 (human-like 42% GC, many homopolymers). The alignment filter removes
those, and deacon removes 2.2M human pairs that `samtools -f 4` (unmapped-only) would have kept.

Outputs per sample:

```
<prefix>.host_depleted.R1.fq.gz
<prefix>.host_depleted.R2.fq.gz
<prefix>.host_depleted.unpaired.fq.gz
<prefix>.deacon.paired.json / .deacon.unpaired.json   # deacon summaries
<prefix>.counts.primary.json / .counts.alignment_filter.json   # samtools read counts
<prefix>.depletion_stats.tsv                          # sample, paired_in/out, unpaired_in/out, frac_kept
```

The script exits non-zero on truncated or corrupt input; it never uploads partial FASTQs.

**Index:** `panhuman-1.k31w15.idx` (HPRC + CHM13v2 + GRCh38, with bacterial and viral minimizers
subtracted; 3.3 GB; ~5 GB RAM). Fetch it with `deacon index fetch panhuman-1` or from
[Zenodo 17288185](https://zenodo.org/records/17288185). It needs deacon ≥ 0.11 (index format v3).

Pinned versions (see `environment.yml`): deacon 0.19.1, samtools 1.25.

## Local use

```bash
mamba env create -f environment.yml
conda activate deacon-host-depletion
scripts/deplete_host.sh -i sample.bam -x panhuman-1.k31w15.idx -p sample -o out -t 16
# CRAM needs the reference it was encoded with:
scripts/deplete_host.sh -i sample.cram -r Homo_sapiens_assembly38.fasta -x panhuman-1.k31w15.idx -p sample -o out
```

Test (synthetic data, no patient reads):

```bash
PATH=$PWD/scripts:$PATH tests/run_synthetic_test.sh
```

## DNAnexus

```bash
dx login
dx select <project>
dnanexus/build.sh --build      # conda-packs samtools+deacon into resources/, then `dx build`
dx upload panhuman-1.k31w15.idx --path /reference/   # once per project/region
dx run deacon_host_depletion -iinput_reads=file-xxx -ideacon_index=file-yyy --destination /host-depleted/
# Batch: dx run deacon_host_depletion --batch-tsv batch.tsv -ideacon_index=file-yyy
```

By default (`-istream=true`, v0.3.0+) the job reads the BAM/CRAM in place through a read-only
[dxfuse](https://github.com/dnanexus/dxfuse) mount (v1.6.1, bundled and checksummed by `build.sh`)
instead of downloading it: reading overlaps with processing and local disk only holds collate temp files
and outputs. Default instance: `azure:mem1_ssd1_x16` (16 cores, 31 GB, 256 GB disk). Samples that keep a
very large fraction of reads (tens of GB of FASTQ) may need a larger disk, e.g. `azure:mem2_ssd1_x16`.

`-istream=false` restores the v0.2 behaviour (download first; the job checks free disk for input × 1.3
and fails fast; use `azure:mem2_ssd1_x16`, or `azure:mem4_ssd1_x32` for BAMs over ~390 GB).

## Cavatica

The image is built, tested and pushed by GitHub Actions (`.github/workflows/image.yml`) to
`ghcr.io/jbk708/deacon-host-depletion`. Pushes to `main` get `:main` and `:sha-xxxxxxx` tags;
a `vX.Y.Z` git tag publishes `:X.Y.Z`. The CWL pins a version tag.

```bash
sbpack cavatica <user>/<project>/deacon-host-depletion cavatica/deacon-host-depletion.cwl
```

Inputs: `input_reads` (CRAM/BAM), `reference` (the Kids First GRCh38 FASTA + `.fai`; required for CRAM),
`deacon_index`, plus optional `threads` (16), `abs_threshold` (2) and `rel_threshold` (0.01).
