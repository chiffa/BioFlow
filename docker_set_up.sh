#!/usr/bin/env bash
# create and set-up home directory
cd /home
mkdir ank
cd ank

# install python build dependencies and system python:
apt-get update
apt-get -y install python
apt-get -y install build-essential
apt-get -y install libsuitesparse-dev
apt-get -y install wget
apt-get -y install git
apt-get -y install curl
apt-get -y install nohup
apt-get -y install lsof

# install and activate Oracle Java 7:
apt-get -y install software-properties-common
apt-get -y install python-software-properties
add-apt-repository ppa:webupd8team/java
apt-get -y update

apt-get -y install oracle-java7-installer
apt-get -y install oracle-java7-set-default

# install neo4j:
cd /home/ank
git clone https://github.com/chiffa/neo4j-community-1.9.6.git
mv neo4j-community-1.9.6 neo4j-yeast

# install mongodb:
cd /home/ank
curl -O https://fastdl.mongodb.org/linux/mongodb-linux-x86_64-3.0.8.tgz
tar -zxvf mongodb-linux-x86_64-3.0.8.tgz
rm mongodb-linux-x86_64-3.0.8.tgz
mv mongodb-linux-x86_64-3.0.8 mongodb
mkdir -p /data/db

# install minicoda
cd /home/ank
wget https://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p /home/ank/miniconda
export PATH="/home/ank/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
rm miniconda.sh

# create and activate test environement:
conda create -q -n test-environement python="2.7" numpy scipy matplotlib
source activate test-environement
conda install python="2.7" cython scikit-learn

# clone the project into the test environement:
mkdir /home/ank/datastore
cd /home/ank
git clone https://github.com/chiffa/BioFlow.git

# install project requirements:
cd /home/ank/BioFlow/
pip install requirements -r requirements.txt

# start-up the databases:
/home/ank/neo4j-yeast/bin/neo4j start
nohup /home/ank/mongodb/bin/mongod &
