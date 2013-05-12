'''
Created on Feb 26, 2013

@author: akucahravy
'''

# TODO: create a simple to use reparser

from TableBuilder import UNIPROT_Prot
from TableBuilder import UNIPROT2GO
from TableBuilder import UNIPROTid2acnum
from TableBuilder import GeneName2UNIPROT
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()

NCBI_TaxID_Interesting=['36329','9606','1773'] #PLAFA taxonomy id is actually 5833

docu=open('/home/akucahravy/Downloads/uniprot_sprot.dat',"r")

i=0
j=0

defpackage={'uniprotid':'err0','full_name':'err1','orf_name':'err2','taxid':'err3','AC_line':' '}
packa=defpackage.copy()
#TODO: redo with a dictionary


################################################
##
## Now, the following part is only to be used on
## the run related to the gene 2 uniprot name table
## filling
##
################################################

Query=session.query(UNIPROT_Prot).all()
    
UNIPROTids=[]
    
for row in Query:
    UNIPROTids.append(str(row.uniprotid))
        
################################################
##
##  Comment it out otherwise
##
################################################

while True:
    i=i+1
    line = docu.readline()
    if not line:
        break
    # print i,line[:-1]
    words=filter(lambda a:a!='', line.split(' '))
    keyword=line[0:2]
    
    if keyword!='//' and keyword!='  ' and len(words)==1:
        bla='bla'
        #print 'key:',keyword,'!!'
        #print i
        #print line
    
    
    if keyword=='//' and packa['taxid'] in NCBI_TaxID_Interesting:
        protein=UNIPROT_Prot(packa['uniprotid'],packa['full_name'],packa['orf_name'],packa['taxid'])
        j=j+1
        #print j,protein
        #session.add(protein)
        for word in filter(lambda a:a!='', packa['AC_line'].split(' '))[1:]:
            sword=word.strip()[:-1]
            Relation2=UNIPROTid2acnum(packa['uniprotid'],sword)
            #print Relation2
            #session.add(Relation2) 
   
        if j%100==0:
            session.commit()
    
    if keyword=='//':
        packa=defpackage.copy()
    
    if keyword=='ID':
        packa['uniprotid']=words[1].strip()
    
    if keyword=='DE':
        if words[1]=='RecName:':
            if 'Full' in words[2]:
                packa['full_name']=line.split('=')[1].split(';')[0].strip()
    #            print packa['full_name']
            else:
                print 'issue here \n',line,'\n','<========================>'
    
    if keyword=='GN':
        if packa['uniprotid'] in UNIPROTids:
            tokens=line[:-2].split(';')
            retokens=[]
            for token in tokens:
                if '=' in token:
                    retokens=retokens+token.split('=')[1].split(', ')
            for retoken in retokens:
                session.add(GeneName2UNIPROT(retoken,packa['uniprotid']))
                print retoken,'|',packa['uniprotid']
        line.split
    
    if keyword=='OX':
        packa['taxid']=line.split('=')[1].split(';')[0].strip()

    if keyword=='DR' and 'GO; GO:' in line and packa['taxid'] in NCBI_TaxID_Interesting:
        Relation=UNIPROT2GO(packa['uniprotid'],line.split(':')[1].split(';')[0])
        #print Relation
        #session.add(Relation)
    
    if keyword=='AC':
        packa['AC_line']=packa['AC_line']+str(line[3:])

session.commit()