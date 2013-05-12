'''
Created on Feb 25, 2013

@author: akucahravy
'''

import TableBuilder
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


lite_engine = create_engine('sqlite:////home/akucahravy/DB/initdb', echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()


docu=open('/home/akucahravy/Downloads/gene_ontology.1_0.obo',"r")

i=0

relations={}

FirstRead=True
LastRead=False

defpackage={'termid':'err0','name':'err1','namespace':'err2','definition':'err3'}
packa=defpackage.copy()

while True:
    i=i+1
    line = docu.readline()
    if not line:
        break
    #print i
        
    attribute = line.split(':')[0].strip()
    
    if attribute == '[Term]':
        if not FirstRead:
            newTerm = TableBuilder.GO_Term(packa['termid'],packa['name'],packa['namespace'],packa['definition'])
            print newTerm
            #TODO: Uncomment this part for Term Loading
            #session.add(newTerm)
        else:
            FirstRead=False
        packa=defpackage.copy()
    
    if attribute =='[Typedef]':
        if not LastRead:
            LastRead=True
            newTerm = TableBuilder.GO_Term(packa['termid'],packa['name'],packa['namespace'],packa['definition'])
            print newTerm
            #TODO: Uncomment this part for Term Loading
            #session.add(newTerm)
    
    if attribute=='id':
        packa['termid']=line.split(':')[-1].strip()
    
    if attribute=='name':
        packa['name']=line.split(':')[-1].strip()
    
    if attribute=='namespace':
        packa['namespace']=line.split(':')[-1].strip()
    
    if attribute=='def':
        packa['definition']=line.split('\"')[1].strip()
            
    if attribute == 'is_a' and not LastRead:
        toGO=line.split(':')[2].split(' ')[0].strip()
        newRel = TableBuilder.GO_Relation(packa['termid'],'is_a',toGO)
        print newRel
        session.add(newRel)
            
    if attribute == 'relationship' and not LastRead:
        relation=line.split(' ')[1]
        
        if relation not in relations.keys():
            relations[relation]=0
        relations[relation]=relations[relation]+1
        toGO=line.split(':')[2].split(' ')[0].strip()
        newRel = TableBuilder.GO_Relation(packa['termid'],relation,toGO)
        print newRel
        session.add(newRel)



# TODO: correct behavior for the last term (action on the block end, not the Term =>open/close behavior)



session.commit()

print relations 
