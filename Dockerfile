FROM mambaorg/micromamba:2.0.5

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app
RUN micromamba run -n base pip install --no-cache-dir ".[ml]"

ENV SPROUT_REF_CACHE=/data/cache \
    SPROUT_REF_BACKBONE=/data/backbone \
    SPROUT_REF_WORK_DIR=/data/runs
EXPOSE 8000
CMD ["sprout-ref", "web", "--host", "0.0.0.0", "--port", "8000"]
