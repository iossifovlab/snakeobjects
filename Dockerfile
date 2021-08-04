FROM bitnami/minideb:buster
MAINTAINER Boris Yamrom <yamrom@cshl.edu>
ADD . /tmp/repo
WORKDIR /tmp/repo
ENV PATH /opt/conda/bin:${PATH}
ENV LANG C.UTF-8
ENV SHELL /bin/bash
RUN install_packages wget curl bzip2 ca-certificates gnupg2 squashfs-tools git
RUN /bin/bash -c "curl -L https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh > miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && \
    rm miniconda.sh"
    
RUN /bin/bash -c "conda install -y -c conda-forge mamba &&                 \
    mamba create -q -y -c iossifovlab                                      \
    -c bioconda -c conda-forge -n snakeobjects snakeobjects --only-deps && \
    source activate snakeobjects &&                                        \
    mamba install -q -y -c conda-forge singularity &&                      \
    conda clean --all -y &&                                                \
    which python &&                                                        \
    pip install .[reports,messaging,google-cloud]"


# Make RUN commands use the new environment:
RUN echo "conda activate snakeobjects" >> ~/.bashrc
RUN echo "export PATH=/workdir/workflow:$PATH" >> ~/.bashrc
ENV PATH /opt/conda/envs/snakeobjects/bin:${PATH}
ENV SO_PIPELINE /workdir/workflow
ENV SO_PROJECT .
ENV SO_CONTAINER "yes"
SHELL ["/bin/bash", "--login", "-c"]




