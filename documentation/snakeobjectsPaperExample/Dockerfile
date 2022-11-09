FROM yamrom/snakeobjects

WORKDIR /app
COPY . /app

RUN /bin/bash -c "conda install -n snakeobjects  \
    -c bioconda -c conda-forge pandas            \
    numpy				         \
    bwakit                                       \
    samtools                                     \
    pysam     &&                                 \
    conda clean --all -y "

# Make RUN commands use the new environment:
RUN echo "conda activate snakeobjects" >> ~/.bashrc
RUN echo "export PATH=/workdir/workflow:$PATH" >> ~/.bashrc
ENV SO_PIPELINE /workdir/workflow
ENV SO_PROJECT .
ENV SO_CONTAINER "yes"
SHELL ["/bin/bash", "--login", "-c"]



