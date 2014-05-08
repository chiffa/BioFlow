"""
funsite here refers to funsite.sdsc.org

This module is no more functional a nd was kept for the legacy reasons
"""
__author__ = 'ank'


# from DBLoader import TableBuilder
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

import PolyPharma.configs as Conf

lite_engine = create_engine(Conf.dbLocation, echo=False)

Session=sessionmaker()
Session.configure(bind=lite_engine)
session=Session()


def PDB2Gene_Parser(GN2PDBname='/home/akucahravy/Downloads/TB-Drugome/pdb_structure_info.csv'):
    '''
    Parses the funsite "pdb_structure_info.***" file.
    Takes as argument the path to the funsite file. 
    Returns a dictionary that maps the Pdb identifiers to the gene names
    '''
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
                
    return PdbR2GN


def drugParser(D2PdbLname='/home/akucahravy/Downloads/TB-Drugome/drug_site_info.csv'):
    '''
    Parses the funsite "drug_site_info.***" file.
    Takes as argument the path to the funsite file. 
    Returns a dictionary that maps drug names to the PDB
        structures they are known to act on.
    
    Fills in the "Drug" SQL table with  the DRUG_Name and LIG_Name
    '''
    
    
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
                session.add(TableBuilder.Drug(code.strip(),drugName))
                LIG2Drug[code.strip()]=drugName
                Drug2LIG[drugName]=code.strip()
            PDBrefs=words[3].split(',')
            D2PdbL[drugName]=[]
            for ref in PDBrefs:
                D2PdbL[drugName].append(ref.strip())
    
    session.commit()
    
    return D2PdbL

def fullPredParse(eHiTsDockName='/home/akucahravy/Downloads/TB-Drugome/SMAP_eHiTs_pdb_structure_info.csv',
                  eHiTsHomoName='/home/akucahravy/Downloads/TB-Drugome/SMAP_eHiTs_homology_model_info.csv'):
    '''
    Fully parses all the files related to prediction of the additional
    docking/homology prediction that points to the relations between the
    drug names and new protein names
    
    Takes as argument the path to the funsite "SMAP_eHiTs_pdb_structure_info.***" 
    and "SMAP_eHiTs_homology_model_info.***" files, in this order.
    
    Fills in the "DRUG2UNIPROT" SQL table with  the LIG_Code and identifier
    '''
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
    #  Ok, so now we will have to zip the structural prediction and the model homology classes
    #  together!
    ##
    
    LIG2Gene=LIG2Dock.copy()                # This and the next line are returned as results, but are also
    LIGGene2Energy=LIFDock2energy.copy()    # But is also used when build to fill in the SQL DRUG2Uniprot table
    
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
    #  Let's check if all the gene names are matching the UNIPROT PROTEINS:
    ##
    
    NotFound=[]
    
    for (LIG,gene) in LIGGene2Energy.keys():
        query = session.query(TableBuilder.GeneName2UNIPROT).\
                            filter(TableBuilder.GeneName2UNIPROT.genename==gene).\
                            all()
        if len(query)>0:
            # we assume tnat the previous reduction step worked fine and 
            # that only one element is present in the return list
            UNIPROTID=query[0].uniprotid
            session.add(TableBuilder.DRUG2UNIPROT(UNIPROTID,LIG,LIGGene2Energy[(LIG,gene)]))
        else:
            NotFound.append((LIG,gene))
    
    session.commit()
    # right now we are just ignoring the not-found category, since it regroups proteins 
    # that are not in SwissProt, but in TrEMBL and thus have very little GO annotation

fullPredParse()