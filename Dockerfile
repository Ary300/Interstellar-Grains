# syntax=docker/dockerfile:1.6
FROM mambaorg/micromamba:1.5.8

USER root
SHELL ["/bin/bash", "-lc"]

WORKDIR /app

COPY environment.yml /app/environment.yml

RUN micromamba create -y -n h2-kmc -f /app/environment.yml && \
    micromamba clean --all --yes

ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV CONDA_DEFAULT_ENV=h2-kmc
ENV PATH=/opt/conda/envs/h2-kmc/bin:$PATH

COPY . /app

CMD ["python", "run_sweep.py"]