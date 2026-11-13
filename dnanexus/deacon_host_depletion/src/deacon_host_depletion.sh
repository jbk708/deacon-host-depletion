#!/bin/bash
# DNAnexus entry point. Tools (samtools, deacon, python) come from the conda-pack tarball
# bundled at resources/opt/deacon-env.tar.gz (built by dnanexus/build.sh), and
# deplete_host.sh is bundled at resources/usr/local/bin/.
set -euo pipefail -x

main() {
    # --- tools ---
    mkdir -p /opt/deacon-env
    tar -xzf /opt/deacon-env.tar.gz -C /opt/deacon-env
    # conda-unpack's shebang is `/usr/bin/env python`, but Ubuntu 24.04 workers only ship python3,
    # so invoke it with the env's own interpreter.
    /opt/deacon-env/bin/python /opt/deacon-env/bin/conda-unpack
    export PATH=/opt/deacon-env/bin:$PATH
    deacon --version; samtools --version | head -1

    threads=$(nproc)
    work=/home/dnanexus/work; mkdir -p "$work/in" "$work/out"

    # --- inputs ---
    in_name=$(dx describe --name "$input_reads")
    # Need room for the input plus collate temp files / outputs; fail fast instead of mid-run.
    in_bytes=$(dx describe "$input_reads" --json | python3 -c 'import json,sys; print(json.load(sys.stdin)["size"])')
    free_bytes=$(df --output=avail -B1 "$work" | tail -1)
    if (( in_bytes * 13 / 10 > free_bytes )); then
        dx-jobutil-report-error "Input is $((in_bytes/1000000000)) GB but only $((free_bytes/1000000000)) GB free; rerun with a larger --instance-type (e.g. azure:mem4_ssd1_x32)." AppError
        exit 1
    fi
    dx download "$input_reads" -o "$work/in/$in_name"
    dx download "$deacon_index" -o "$work/in/deacon.idx"
    ref_opt=()
    if [[ -n "${reference:-}" ]]; then
        ref_name=$(dx describe --name "$reference")
        dx download "$reference" -o "$work/in/$ref_name"
        ref_opt=(-r "$work/in/$ref_name")
    fi
    prefix=${output_prefix:-$(echo "$in_name" | sed -E 's/\.(bam|cram)$//')}

    df -h "$work"

    # --- run ---
    deplete_host.sh -i "$work/in/$in_name" -x "$work/in/deacon.idx" -p "$prefix" \
        -o "$work/out" -t "$threads" -a "$abs_threshold" -R "$rel_threshold" \
        -T "$work/tmp" "${ref_opt[@]}"

    # free the input before uploading
    rm -f "$work/in/$in_name"

    # --- outputs ---
    r1=$(dx upload "$work/out/$prefix.host_depleted.R1.fq.gz" --brief)
    r2=$(dx upload "$work/out/$prefix.host_depleted.R2.fq.gz" --brief)
    up=$(dx upload "$work/out/$prefix.host_depleted.unpaired.fq.gz" --brief)
    st=$(dx upload "$work/out/$prefix.depletion_stats.tsv" --brief)
    dx-jobutil-add-output r1 "$r1" --class=file
    dx-jobutil-add-output r2 "$r2" --class=file
    dx-jobutil-add-output unpaired "$up" --class=file
    dx-jobutil-add-output stats "$st" --class=file
    for j in "$work"/out/*.deacon.*.json; do
        dx-jobutil-add-output deacon_summaries "$(dx upload "$j" --brief)" --class=array:file
    done
}
