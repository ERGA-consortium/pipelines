FROM continuumio/miniconda3
ADD environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
# Pull the environment name out of the environment.yml
RUN echo "source activate gaas" > ~/.bashrc
ENV PATH /opt/conda/envs/protexcluder/bin:$PATH
