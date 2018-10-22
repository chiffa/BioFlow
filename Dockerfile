FROM ubuntu:18.04
MAINTAINER "Andrei Kucharavy <ank@andreikucharavy.com>"
# change to continuumio/anaconda

RUN export DEBIAN_FRONTEND=noninteractive
RUN echo 'debconf debconf/frontend select Noninteractive' | debconf-set-selections

# create and set-up home directory:
RUN cd /home
RUN mkdir /home/ank
RUN cd /home/ank

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
RUN cd /home/ank
ADD https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh miniconda.sh
RUN bash miniconda.sh -b -p /home/ank/miniconda
ENV PATH="/home/ank/miniconda/bin:${PATH}"
RUN hash -r
RUN conda config --set always_yes yes --set changeps1 no
RUN conda update -q conda
RUN rm miniconda.sh

# create and activate conda environment
RUN conda create -q -n run-environement python="2.7" numpy=1.9 scipy=0.19 matplotlib=1.4
RUN /bin/bash -c "source activate run-environement"
RUN conda install python="2.7" cython=0.22 scikit-learn=0.16

# clone the project into the test environement:
ADD https://github.com/chiffa/BioFlow/archive/master.zip BioFlow.zip
RUN unzip BioFlow.zip

# install project requirements:
RUN cd /BioFlow-master/; pip install requirements -r requirements.txt


