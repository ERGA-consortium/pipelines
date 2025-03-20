FROM mambaorg/micromamba
ADD environment.yml /tmp/environment.yml
RUN micromamba create -f /tmp/environment.yml -y && micromamba clean --all --yes
RUN echo "source activate biopython-aiohttp" > ~/.bashrc
ENV PATH /opt/conda/envs/biopython-aiohttp/bin:$PATH
