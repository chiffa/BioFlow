'''
Created on Mar 14, 2013

@author: akucahravy
'''



##
# SMAP ehits file methionned here can be downloaded from funsite.sdsc.edu
##

from DBLoader import TableBuilder
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()


GN2PDBname='/home/akucahravy/Downloads/TB-Drugome/pdb_structure_info.csv'
GN2PDB=open(GN2PDBname,'r')

PdbR2GN={}

i=0
while True:
    i=i+1
    line = GN2PDB.readline()
    if not line:
        break
    
    if i>2:
        words=line.split('\t')
        geneName=words[0].strip()
        PDBrefs=words[3].replace(', "','').replace('"', '').strip().split(',')
        for ref in PDBrefs:
            PdbR2GN[ref.strip()]=geneName


###
#
#  Parse the drug sheet to recover the drugNames and finish building the DB structure
#
###

D2PdbLname='/home/akucahravy/Downloads/TB-Drugome/drug_site_info.csv'
D2PdbLfile=open(D2PdbLname,'r')

D2PdbL={}
LIG2Drug={}
Drug2LIG={}

i=0

while True:
    i=i+1
    line = D2PdbLfile.readline()
    if not line:
        break

    if i>2:
        words=line.split('\t')
        drugName=words[0].strip()
        ligandcode=words[1].split(',')
        for code in ligandcode:
            LIG2Drug[code.strip()]=drugName
            Drug2LIG[drugName]=code.strip()
        PDBrefs=words[3].split(',')
        D2PdbL[drugName]=[]
        for ref in PDBrefs:
            D2PdbL[drugName].append(ref.strip())


##
# Now we can do drug- > LIG Key -> PDB id -> Gene name -> UNIPROT_ID
# and then root from Uniprot ID directly to the annotation
# space and to the BioPax format
##


######
#
#  This part is unnecessary, because all the known drug targets have been used to position
#  it against it's original disease and are not necessary present in the organisms
#  we are interested in.
#
######
#LIG2Known_Gene={}
#for Drug in D2PdbL.keys():
#    print '<============>'
#    print Drug
#    LIG=Drug2LIG[Drug]
#    print LIG
#    Genes=[]
#    for PDBid in D2PdbL[Drug]:
#        print PDBid
#        Gene=PdbR2GN[PDBid]
#        print Gene
#        Genes.append(Gene)
#    LIG2Known_Gene[LIG]=Genes[:]


eHiTsDockName='/home/akucahravy/Downloads/TB-Drugome/SMAP_eHiTs_pdb_structure_info.csv'
eHiTsDockFile=open(eHiTsDockName,'r')

LIG2Dock={}
LIFDock2energy={}

i=0

while True:
    i=i+1
    line = eHiTsDockFile.readline()
    if not line:
        break

    if i>2:
        words=line.split('\t')
        #targetPDB=words[0] # if you want the PDB identifier rather then gene
        Gene=words[2].replace('"', '').strip().split(', ')
        DrugLIG=words[3].strip()
        #Ehits score is assumed to be proportional to the energy.
        EhitsEnergy=words[8].strip()
        # we set here a default Ehits Score. the value chosen is the average of lowest quartile and minimum: -4.41 
        if EhitsEnergy=='x':
            EhitsEnergy='-4.41'
        EhitsEnergy=float(EhitsEnergy)
        if DrugLIG not in LIG2Dock.keys():
            LIG2Dock[DrugLIG]=[]
        LIG2Dock[DrugLIG]=LIG2Dock[DrugLIG]+Gene
        for gene in Gene:
            if (DrugLIG,gene) not in LIFDock2energy.keys():
                LIFDock2energy[(DrugLIG,gene)]=[]
            LIFDock2energy[(DrugLIG,gene)].append(EhitsEnergy)
        
# this part performs the redundancy reductuion in the genes attached to a drug through
# different models and calculates the average eHiTs score

for LIG in LIG2Dock.keys():
    LIG2Dock[LIG]=list(set(LIG2Dock[LIG]))
    for gene in LIG2Dock[LIG]:
        LIFDock2energy[(LIG,gene)]=sum(LIFDock2energy[(LIG,gene)])/float(len(LIFDock2energy[(LIG,gene)]))

eHiTsHomoName='/home/akucahravy/Downloads/TB-Drugome/SMAP_eHiTs_homology_model_info.csv'
eHiTsHomoFile=open(eHiTsHomoName,'r')

LIG2Homo={}
LIFHomo2energy={}

i=0

while True:
    i=i+1
    line = eHiTsHomoFile.readline()
    if not line:
        break

    if i>2:
        words=line.split('\t')
        #targetPDB=words[0] # if you want the PDB identifier rather then gene
        Gene=words[1].replace('"', '').strip().split(', ')
        DrugLIG=words[3].strip()
        #Ehits score is assumed to be proportional to the energy.
        EhitsEnergy=words[5].strip()
        # we set here a default Ehits Score. the value chosen is the average of lowest quartile and minimum: -4.41 
        if EhitsEnergy=='x':
            EhitsEnergy='-4.41'
        EhitsEnergy=float(EhitsEnergy)
        if DrugLIG not in LIG2Homo.keys():
            LIG2Homo[DrugLIG]=[]
        LIG2Homo[DrugLIG]=LIG2Homo[DrugLIG]+Gene
        for gene in Gene:
            if (DrugLIG,gene) not in LIFHomo2energy.keys():
                LIFHomo2energy[(DrugLIG,gene)]=[]
            LIFHomo2energy[(DrugLIG,gene)].append(EhitsEnergy)
        
# this part performs the redundancy reductuion in the genes attached to a drug through
# different models and calculates the average eHiTs score

for LIG in LIG2Homo.keys():
    LIG2Homo[LIG]=list(set(LIG2Homo[LIG]))
    for gene in LIG2Homo[LIG]:
        LIFHomo2energy[(LIG,gene)]=sum(LIFHomo2energy[(LIG,gene)])/float(len(LIFHomo2energy[(LIG,gene)]))            
        
##
#
#  Ok, so now we will have to zip the structural prediction and the model homology classes
#  together!
#
##


LIG2Gene=LIG2Dock.copy() #This and the next line have to become system-wide variable
LIGGene2Energy=LIFDock2energy.copy()

for LIG in LIG2Homo.keys():
    if LIG not in LIG2Gene.keys():
        LIG2Gene[LIG]=[]
    LIG2Gene[LIG]=LIG2Gene[LIG]+LIG2Homo[LIG]
    LIG2Gene[LIG]=list(set(LIG2Gene[LIG]))
    for gene in LIG2Homo[LIG]:
        LIGGene2Energy[(LIG,gene)]=LIFHomo2energy[(LIG,gene)]
        
#for LIG in LIG2Gene.keys():
#    print LIG,LIG2Gene[LIG]
#    for gene in LIG2Gene[LIG]:
#        print (LIG, gene), LIGGene2Energy[(LIG, gene)]
        

##
#
#  Let's check if all the gene names are matching the UNIPROT PROTEINS:
#
##

NotFound=[]

i=0
for values in LIG2Gene.values():
    for val in values:
        i+=1
        query = session.query(TableBuilder.GeneName2UNIPROT).\
                    filter(TableBuilder.GeneName2UNIPROT.genename==val).\
                    all()
        if len(query)>0:
            print query
        else:
            NotFound.append(val)
print 'Not found : ', len(NotFound), 'out of ', len(NotFound)+i
for value in NotFound:
    print value
    



