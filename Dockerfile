FROM bitnami/minideb:bullseye

RUN install_packages wget curl bzip2 ca-certificates gnupg2 squashfs-tools git

RUN /bin/bash -c "curl -L https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh"
ENV PATH /opt/conda/bin:${PATH}

RUN /bin/bash -c "conda install -y -c conda-forge mamba                            && \
    mamba create -q -y -n snakeobjects -c iossifovlab -c bioconda -c conda-forge snakeobjects singularity && \
    conda clean --all -y                                                           && \
    which python"

ENV PATH /opt/conda/envs/snakeobjects/bin:${PATH}

