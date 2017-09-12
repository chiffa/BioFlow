FROM ubuntu:14.04

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

# install and activate Oracle Java 7:
RUN apt-get -yq install software-properties-common
RUN apt-get -yq install python-software-properties

# install and activate Oracle Java 7:
RUN echo "deb http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" | tee -a /etc/apt/sources.list
RUN echo "deb-src http://ppa.launchpad.net/webupd8team/java/ubuntu trusty main" | tee -a /etc/apt/sources.list
RUN echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys EEA14886 && apt-get update

RUN wget https://www.dropbox.com/sh/gxeklrzkq58ydsf/AADF4dimrwsmsprUxgm2Iwn8a/jdk-7u80-linux-x64.tar.gz?dl=1
RUN mv jdk-7u80-linux-x64.tar.gz?dl=1 jdk-7u80-linux-x64.tar.gz
RUN mkdir /var/cache/oracle-jdk7-installer/
RUN sudo mv jdk-7u80-linux-x64.tar.gz /var/cache/oracle-jdk7-installer/

RUN apt-get install -y curl dnsutils oracle-java7-installer ca-certificates
RUN apt-get -yq install oracle-java7-set-default

# install neo4j:
RUN cd /home/ank
RUN git clone https://github.com/chiffa/neo4j-community-1.9.6.git
RUN mv neo4j-community-1.9.6 neo4j-yeast

# install mongodb:
RUN cd /home/ank
RUN curl -O https://fastdl.mongodb.org/linux/mongodb-linux-x86_64-3.0.8.tgz
RUN tar -zxvf mongodb-linux-x86_64-3.0.8.tgz
RUN rm mongodb-linux-x86_64-3.0.8.tgz
RUN mv mongodb-linux-x86_64-3.0.8 mongodb
RUN mkdir -p /data/db

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


# you will need to connect to the container and run those commands to get the databases up
# /neo4j-yeast/bin/neo4j start
# /mongodb/bin/mongod &


