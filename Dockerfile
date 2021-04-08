FROM ubuntu:20.04
MAINTAINER "Andrei Kucharavy <ank@andreikucharavy.com>"
# change to continuumio/anaconda

RUN export DEBIAN_FRONTEND=noninteractive
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

RUN cd /home

# install python build dependencies and system python:
RUN apt-get update
RUN apt-get -yq install python
RUN apt-get -yq install build-essential
RUN apt-get -yq install libsuitesparse-dev
RUN apt-get -yq install wget
RUN apt-get -yq install git
RUN apt-get -yq install curl
RUN apt-get -yq install unzip
RUN apt-get -yq install lsof
RUN apt-get update
RUN apt-get -yq install libsm6 libxrender1 libfontconfig1 libglib2.0-0

# install minicoda
ADD https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh miniconda.sh
RUN bash miniconda.sh -b -p /miniconda
ENV PATH="/miniconda/bin:${PATH}"
RUN hash -r
RUN conda config --set always_yes yes --set changeps1 no
RUN conda update -q conda
RUN rm miniconda.sh

# create and activate conda environment
RUN conda create -q -n run-environement python="3.7" numpy scipy matplotlib
RUN /bin/bash -c "source activate run-environement"
RUN conda install python="3.7" cython scikit-learn

# clone the project into the test environement:
ADD https://github.com/chiffa/BioFlow/archive/master.zip BioFlow.zip
RUN unzip BioFlow.zip
RUN rm BioFlow.zip
RUN apt-get install -yq nano

# install project requirements:
RUN which pip
RUN cat /BioFlow-master/requirements.txt
RUN cd /BioFlow-master/; pip install -r requirements.txt
