'''
Created on Apr 9, 2013

@author: akucahravy
'''

from DBLoader import TableBuilder as TB
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


lite_engine = create_engine('sqlite:////home/akucahravy/DB/initdb', echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()


infile='/home/akucahravy/Downloads/pdbtosp.txt'
outfile='/home/akucahravy/pdb2UNIPROT_reinc.txt'

inf=open(infile,'r')

selline='Swiss-Prot entry name(s)'

dico={}

k=0
currWord=''
while True:
    line=inf.readline()
    if not line:
        break
    k=k+1
    if k>25:
        words=line.strip().split(' ')
        words=filter(lambda a:a !='',words)
        if 'X-ray' in words:
            currWord=words[0]
            dico[currWord]=[]
            for i in range(0,len(words[4:])/2):
                dico[currWord].append(words[4+2*i])
        else :
            if 'NMR' in words:
                currWord=words[0]
                dico[currWord]=[]
                for i in range(0,len(words[3:])/2):
                    dico[currWord].append(words[3+2*i])
            else:
                for i in range(0,len(words)/2):
                    dico[currWord].append(words[2*i])
   
for key,llist in dico.items():
    for val in llist:
        session.add(TB.PDB2UNIPROT(key,val))

session.commit()