FROM snakemake/snakemake

WORKDIR /app
COPY . /app

#RUN /bin/bash -c "conda install -n snakemake -c iossifovlab \
#    -c bioconda -c conda-forge snakeobjects && \
#    conda clean --all -y "

RUN /bin/bash -c "cd /app && pip install . && \
    conda clean --all -y "

# Make RUN commands use the new environment:
RUN echo "conda activate snakemake" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]





