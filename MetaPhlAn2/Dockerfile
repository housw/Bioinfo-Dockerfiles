FROM shengwei/ubuntu1604py27:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# install metaphlan2
RUN conda install -c bioconda numpy matplotlib scipy metaphlan2 && \
    conda clean -ya 

#  'axisbg' in 'artist.py' has changed to 'facecolor' from matplotlib v2.2.3
RUN /opt/conda/bin/pip uninstall -y matplotlib && /opt/conda/bin/pip install matplotlib==2.1.0

#  Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

