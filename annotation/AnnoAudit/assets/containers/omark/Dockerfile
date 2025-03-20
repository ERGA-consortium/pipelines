FROM continuumio/miniconda3
ADD environment.yml /tmp/environment.yml
RUN conda install --yes nomkl
RUN conda env create -f /tmp/environment.yml && conda clean -afy
RUN echo "source activate omark" > ~/.bashrc
ENV PATH /opt/conda/envs/omark/bin:$PATH
