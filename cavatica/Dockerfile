# Image for the Cavatica (Seven Bridges) CWL tool.
# Built and pushed by .github/workflows/image.yml to ghcr.io/jbk708/deacon-host-depletion.
# Local build (from the repo root): docker build -f cavatica/Dockerfile -t deacon-host-depletion .
FROM mambaorg/micromamba:1.5.10

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

USER root
COPY scripts/deplete_host.sh /usr/local/bin/deplete_host.sh
RUN chmod 755 /usr/local/bin/deplete_host.sh
ENV PATH=/opt/conda/bin:$PATH
# Seven Bridges runs commands with `sh -c`; reset the micromamba entrypoint so PATH above is used directly.
ENTRYPOINT []
