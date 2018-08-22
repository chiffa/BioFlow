FROM ubuntu:14.04
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
#RUN apt-get -yq install nohup
RUN apt-get -yq install lsof
RUN apt-get -yq install libsm6 libxrender1 libfontconfig1 libglib2.0-0

# install minicoda
RUN cd /home/ank
RUN wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -p /home/ank/miniconda
ENV PATH="/home/ank/miniconda/bin:${PATH}"
RUN hash -r
RUN conda config --set always_yes yes --set changeps1 no
RUN conda update -q conda
RUN rm miniconda.sh

# create and activate test environement:
RUN conda create -q -n test-environement python="2.7" numpy scipy matplotlib
RUN /bin/bash -c "source activate test-environement"
RUN conda install python="2.7" cython scikit-learn

# clone the project into the test environement:
RUN mkdir /home/ank/datastore
RUN cd /home/ank
RUN git clone https://github.com/chiffa/BioFlow.git

# install project requirements:
RUN /bin/bash -c "cd /BioFlow/; pip install requirements -r requirements.txt"

# TODO: is this still true?
# you will need to connect to the container and run those commands to get the databases up
# /neo4j-yeast/bin/neo4j start
# /mongodb/bin/mongod &


