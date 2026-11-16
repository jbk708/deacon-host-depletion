#!/usr/bin/env bash
# Stage resources for the DNAnexus applet, then (optionally) build it.
#   dnanexus/build.sh            # stage resources only
#   dnanexus/build.sh --build    # stage + `dx build` into the current dx project
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
repo=$(dirname "$here")
app=$here/deacon_host_depletion
res=$app/resources
MAMBA=${MAMBA:-mamba}

mkdir -p "$res/opt" "$res/usr/local/bin"
install -m 755 "$repo/scripts/deplete_host.sh" "$res/usr/local/bin/deplete_host.sh"

envdir=$(mktemp -d)/env
"$MAMBA" env create -y -p "$envdir" -f "$repo/environment.yml"
"$MAMBA" install -y -p "$envdir" -c conda-forge conda-pack
"$envdir/bin/conda-pack" -p "$envdir" -o "$res/opt/deacon-env.tar.gz" --force \
  --exclude 'lib/python*/site-packages/conda_pack*'
chmod 644 "$res/opt/deacon-env.tar.gz"
ls -lh "$res/opt/deacon-env.tar.gz"

if [[ ${1:-} == --build ]]; then
  dx build "$app" --overwrite --destination "${DX_APPLET_DEST:-/applets/}"
fi
