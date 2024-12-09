# Install neo4j

sudo apt install openjdk-17-jre-headless
wget -O - https://debian.neo4j.com/neotechnology.gpg.key | sudo gpg --dearmor -o /etc/apt/keyrings/neotechnology.gpg
echo 'deb [signed-by=/etc/apt/keyrings/neotechnology.gpg] https://debian.neo4j.com stable latest' | sudo tee -a /etc/apt/sources.list.d/neo4j.list
sudo apt-get update
apt list -a neo4j
sudo apt-get install neo4j=1:5.25.1
neo4j-admin dbms set-initial-password <NEOPASS>
sudo systemctl start neo4j
sudo systemctl status neo4j


#Install mongodb

sudo apt-get install gnupg curl
curl -fsSL https://www.mongodb.org/static/pgp/server-8.0.asc |    sudo gpg -o /usr/share/keyrings/mongodb-server-8.0.gpg    --dearmor
lsb_release -a
echo "deb [ arch=amd64,arm64 signed-by=/usr/share/keyrings/mongodb-server-8.0.gpg ] https://repo.mongodb.org/apt/ubuntu jammy/mongodb-org/8.0 multiverse" | sudo tee /etc/apt/sources.list.d/mongodb-org-8.0.list
sudo apt-get update
sudo apt-get install -y mongodb-org
sudo systemctl start mongod
sudo systemctl status mongod



# Install Miniconda:

wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" -O ~/miniforge.sh
bash ~/miniforge.sh -b -p ~/miniforge3
source ~/miniforge3/etc/profile.d/conda.sh
conda activate
conda init
exit

<re-log into the terminal>


# Intall BioFlow

git clone https://github.com/chiffa/BioFlow.git
cd BioFlow/
sudo apt-get -y install libsuitesparse-dev
conda create -q -n test-environment python=3.11 numpy scipy matplotlib
conda activate bioflowenv
conda create -q -n bioflowenv python=3.11 numpy scipy matplotlib
conda install python=3.11 cython scikit-learn coverage pylint
pip install coveralls pyflakes pep8-naming mccabe flake8 codecov pyyaml
conda install -c conda-forge scikit-sparse
pip install -r requirements.txt

# run the test analysis pipeline

cd BioFlow/
conda activate bioflowenv
export NEOPASS=<NEOPASS>
python -m bioflow.analysis_pipeline_example