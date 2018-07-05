# shengwei/datascience-notebook 
#
# VERSION       latest

FROM jupyter/datascience-notebook
MAINTAINER Shengwei Hou, housw2010@gmail.com


USER $NB_UID

# update root env 
COPY packages/packages.yml .
RUN conda env update --name root --file packages.yml
RUN rm *.yml

# install julia packages
RUN julia -e 'Pkg.init()' && \
    julia -e 'Pkg.update()' && \
    julia -e 'Pkg.add("FastaIO")' && \
    julia -e 'Pkg.add("ArgParse")'
