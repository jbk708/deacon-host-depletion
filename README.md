# deacon-host-depletion

Host (human) read depletion for aligned WGS BAM/CRAM files using
[deacon](https://github.com/bede/deacon), packaged for two platforms:

| Platform | Data | Wrapper |
|---|---|---|
| DNAnexus (St. Jude Cloud) | St. Jude WGS BAMs (`azure:westus`) | `dnanexus/deacon_host_depletion/` applet |
| Cavatica (Seven Bridges) | CBTN / PBTA WGS CRAMs | `cavatica/deacon-host-depletion.cwl` + Docker image |

Both call the same script, `scripts/deplete_host.sh`.

## What it does

1. Streams **all** reads (primary alignments, duplicates and QC-fail reads included) out of the
   BAM/CRAM with `samtools collate -f | samtools fastq`, so depletion does not rely on the
   original alignment.
2. Pipes interleaved pairs into `deacon filter --deplete --interleaved --check-pairs`, which writes
   R1/R2 files with mates kept in sync.
3. Runs singletons (mate missing) and unpaired reads through deacon separately.

Outputs per sample:

```
<prefix>.host_depleted.R1.fq.gz
<prefix>.host_depleted.R2.fq.gz
<prefix>.host_depleted.unpaired.fq.gz
<prefix>.deacon.paired.json / .deacon.unpaired.json   # deacon summaries
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

Default instance: `azure:mem1_ssd1_x16` (`azure:westus`) or `mem1_ssd1_v2_x16` (`aws:us-east-1`).
Confirm the local disk is larger than the biggest BAM plus collate temp space before running at scale.

## Cavatica

The image is built, tested and pushed by GitHub Actions (`.github/workflows/image.yml`) to
`ghcr.io/jbk708/deacon-host-depletion`. Pushes to `main` get `:main` and `:sha-xxxxxxx` tags;
a `vX.Y.Z` git tag publishes `:X.Y.Z`. The CWL pins a version tag.

```bash
sbpack cavatica <user>/<project>/deacon-host-depletion cavatica/deacon-host-depletion.cwl
```

Inputs: `input_reads` (CRAM/BAM), `reference` (the Kids First GRCh38 FASTA + `.fai`; required for CRAM),
`deacon_index`, plus optional `threads` (16), `abs_threshold` (2) and `rel_threshold` (0.01).
