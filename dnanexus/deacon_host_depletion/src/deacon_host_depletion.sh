#!/bin/bash
# DNAnexus entry point. Tools (samtools, deacon, python) come from the conda-pack tarball
# bundled at resources/opt/deacon-env.tar.gz, and deplete_host.sh and dxfuse are bundled at
# resources/usr/local/bin/ (all staged by dnanexus/build.sh).
#
# stream=true (default): the input BAM/CRAM is read in place through a read-only dxfuse mount, so the
# job never copies it to local disk; reading overlaps with processing, and the instance only needs disk
# for samtools collate temp files and the outputs.
# stream=false: download the input first (the v0.2 behaviour; needs disk for input x 1.3).
set -euo pipefail -x

# Locate a file inside a dxfuse mount by its file ID (dxfuse exposes it as the `base.id` xattr). Paths
# are not guessed because dxfuse renames files with duplicate names or slashes.
find_in_mount() {  # mount_root folder name file_id
    python3 - "$@" <<'PY'
import os, sys
root, folder, name, fid = sys.argv[1:]
for proj in sorted(os.listdir(root)):
    top = os.path.join(root, proj, folder.lstrip("/"))
    for dirpath, dirnames, files in os.walk(top):
        if os.path.relpath(dirpath, top).count(os.sep) >= 1:
            dirnames[:] = []  # duplicate-name files sit at most two levels below the folder
        if name in files:
            p = os.path.join(dirpath, name)
            try:
                if os.getxattr(p, "base.id").decode() == fid:
                    print(p)
                    sys.exit(0)
            except OSError:
                pass
sys.exit(f"{fid} ({folder}/{name}) not found under {root}")
PY
}

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
    in_json=$(dx describe "$input_reads" --json)
    in_name=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["name"])' <<<"$in_json")
    in_id=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$in_json")
    in_proj=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["project"])' <<<"$in_json")
    in_folder=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["folder"])' <<<"$in_json")
    in_bytes=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["size"])' <<<"$in_json")

    if [[ "${stream:-true}" == "true" ]]; then
        mnt=/home/dnanexus/dxmnt; mkdir -p "$mnt"
        sudo_cmd=(); [[ $(id -u) -ne 0 ]] && sudo_cmd=(sudo -E)
        "${sudo_cmd[@]}" chmod u+rw /dev/fuse || true
        "${sudo_cmd[@]}" dxfuse -readOnly "$mnt" "$in_proj"
        in_path=$(find_in_mount "$mnt" "$in_folder" "$in_name" "$in_id")
        echo "streaming $in_path ($((in_bytes/1000000000)) GB) via dxfuse"
    else
        # Need room for the input plus collate temp files / outputs; fail fast instead of mid-run.
        free_bytes=$(df --output=avail -B1 "$work" | tail -1)
        if (( in_bytes * 13 / 10 > free_bytes )); then
            dx-jobutil-report-error "Input is $((in_bytes/1000000000)) GB but only $((free_bytes/1000000000)) GB free; rerun with a larger --instance-type (e.g. azure:mem4_ssd1_x32) or stream=true." AppError
            exit 1
        fi
        dx download "$input_reads" -o "$work/in/$in_name"
        in_path="$work/in/$in_name"
    fi
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
    deplete_host.sh -i "$in_path" -x "$work/in/deacon.idx" -p "$prefix" \
        -o "$work/out" -t "$threads" -a "$abs_threshold" -R "$rel_threshold" \
        -T "$work/tmp" "${ref_opt[@]}"

    # free the input / mount before uploading
    if [[ "${stream:-true}" == "true" ]]; then
        "${sudo_cmd[@]}" fusermount -u "$mnt" || true
    else
        rm -f "$work/in/$in_name"
    fi

    # --- outputs ---
    r1=$(dx upload "$work/out/$prefix.host_depleted.R1.fq.gz" --brief)
    r2=$(dx upload "$work/out/$prefix.host_depleted.R2.fq.gz" --brief)
    up=$(dx upload "$work/out/$prefix.host_depleted.unpaired.fq.gz" --brief)
    st=$(dx upload "$work/out/$prefix.depletion_stats.tsv" --brief)
    dx-jobutil-add-output r1 "$r1" --class=file
    dx-jobutil-add-output r2 "$r2" --class=file
    dx-jobutil-add-output unpaired "$up" --class=file
    dx-jobutil-add-output stats "$st" --class=file
    for j in "$work"/out/*.json; do
        dx-jobutil-add-output deacon_summaries "$(dx upload "$j" --brief)" --class=array:file
    done
}
