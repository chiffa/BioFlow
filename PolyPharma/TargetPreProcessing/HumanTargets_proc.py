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


infile='/home/akucahravy/Downloads/interaction.dat'
outfile='/home/akucahravy/Downloads/interaction_corr.dat'

secEffFile=open('/home/akucahravy/Downloads/meddra_adverse_effects.tsv','r')

inf=open(infile,'r')
outf=open(outfile,'w')

finlist=['sulindac', 'cytosine arabinoside', 'mupirocin',
        'tamoxifen', 'daunorubicin', 'triiodothyronine',
        'saquinavir', 'delavirdine', 'dopamine', 'ampicillin',
        'rimantadine', 'retinoic acid', 'isoflurane', 'trimethoprim',
        'vardenafil', 'felodipine', 'tobramycin', 'rapamycin',
        'norethisterone', 'celecoxib', 'tacrine', 'paromomycin',
        'montelukast', 'raltitrexed', 'paclitaxel', 'progesterone',
        'raloxifene', 'prochlorperazine', 'sunitinib', 'amikacin',
        'furosemide', 'nevirapine', 'diclofenac', 'imatinib',
        'ibuprofen', 'methotrexate', 'diflunisal', 'trifluoperazine',
        'nilotinib', 'mefenamic acid', 'caffeine', 'doxepin',
        'estradiol', 'fluconazole', 'bicalutamide', 'cladribine',
        'thyroxine', 'clofarabine', 'halothane', 'ethacrynic acid',
        'iodipamide', 'nelfinavir', 'mifepristone', "2'-deoxycoformycin",
        'pemetrexed', 'dexamethasone', 'miconazole', 'testosterone', 
        'posaconazole', 'bimatoprost', 'gefitinib', 'efavirenz', 
        'ibandronate', 'trimetrexate', 'timolol', 'oxytetracycline', 
        'doxycycline', 'dobutamine', 'rosiglitazone', 'erythromycin', 
        'leucovorin', 'colchicine', 'vinorelbine', 'econazole', 'chloroquine', 
        'pioglitazone', 'ritonavir', 'chloramphenicol', 'lovastatin', 
        'salbutamol', 'penciclovir', 'indinavir', 'docetaxel', 'rifampicin', 
        'sildenafil', 'spironolactone', 'adenosine', 'mycophenolic acid', 
        'dasatinib', 'pyrimethamine', 'aliskiren', 'tadalafil', 'erlotinib', 
        'diazepam', 'acyclovir', 'd-penicillamine', 'digoxin', 'ketoconazole', 
        'atazanavir', 'amantadine', 'tetrahydrobiopterin', 'amprenavir', 
        'flurbiprofen', 'tetracycline', 'fludrocortisone', 'amphetamine', 
        'propofol', 'isoproterenol', 'bexarotene', 'darunavir', 
        'mitoxantrone', 'tacrolimus', 'sorafenib', 'zidovudine', 
        'salicylic acid', 'levonorgestrel', 'indomethacin']

mapFile={'mycophenolic':'mycophenolic acid',
         'pentostatin':"2'-deoxycoformycin",
         'cytarabine':'cytosine arabinoside',
         'ethacrynic':'ethacrynic acid',
         'norethindrone':'norethisterone',
         'sirolimus':'rapamycin',
         'mefenamic':'mefenamic acid',
         'levothyroxine':'thyroxine',
         'penicillamine':'d-penicillamine',
         'rifampin':'rifampicin',
         'alitretinoin':'retinoic acid',
         'liothyronine':'triiodothyronine',
         'aciclovir':'acyclovir',
         'db00586':'diclofenac',
         'norgestrel':'levonorgestrel',
         'salicyclic':'salicylic acid',
         'vinblastine':'vinorelbine',
         'kanamycin':'amikacin',
         }

revMapFile = {v:k for k,v in mapFile.items()}

def report_uniprot(UNIP_ID):
    result=[]
    accq=session.query(TB.UNIPROT_Prot).filter(TB.UNIPROT_Prot.uniprotid==UNIP_ID).one()
    result.append(UNIP_ID)
    result.append(str(accq.full_name))
    geneq=session.query(TB.GeneName2UNIPROT).filter(TB.GeneName2UNIPROT.uniprotid==UNIP_ID).all()
    buffer1=''
    for row in geneq:
        buffer1+=str(row.genename)+' '
    result.append(buffer1)
    accnumq=session.query(TB.UNIPROTid2acnum).filter(TB.UNIPROTid2acnum.uniprotid==UNIP_ID).all()
    buffer2=''
    for row in accnumq:
        buffer2+=str(row.acnum)+' '
    result.append(buffer2)
    
    return result
        
    
def Li_report_buider():
    i=0
    while True:
        line=inf.readline()
        if not line:
            break
        i+=1
        
        words=line.strip().split('\t')
        quer=session.query(TB.PDB2UNIPROT).\
                            filter(TB.PDB2UNIPROT.pdbid==words[0]).\
                            all()
        UNIPROTS=[]
        for row in quer:
            UNIPROTS.append(str(row.uniprotid))
            
        Buffer=''
        for word in words:
            Buffer+=word+'\t'
        
        for uniprot in UNIPROTS:
            count=session.query(TB.UNIPROT_Prot).filter(TB.UNIPROT_Prot.uniprotid==uniprot).count()
            if count>0:
                lllist=report_uniprot(uniprot)
                Buffer+='\n\t\t\t'
                for elt in lllist:
                    Buffer+=elt+'\t'
        
        outf.write(Buffer)
        print Buffer

drugname=[]
                        
def Human_proteome_Loader():
    i=0
    while True:
        line=inf.readline()
        if not line:
            break
        i+=1
        
        words=line.strip().split('\t')
        quer=session.query(TB.PDB2UNIPROT).\
                            filter(TB.PDB2UNIPROT.pdbid==words[0]).\
                            all()
        
        for row in quer:
            if '_HUMAN' in row.uniprotid:
                if words[1].lower() in mapFile.keys():
                    drugname.append(mapFile[words[1].lower()])
                else:
                    drugname.append(words[1].lower())
            session.add(TB.DRUG2HUMUNIPROT(row.uniprotid,words[1].lower(),words[2]))

SecEffDrugNames=[]
SecEffIDs=[]

def SecEff_Drugs_Scan():
    i=0
    while True:
        line=secEffFile.readline()
        if not line:
            break
        i+=1
        
        words=line.strip().split('\t')
        SecEffDrugNames.append(words[3].lower())
        SecEffIDs.append((words[6],words[7]))

        if words[5]=='PT' and words[3].lower() in finlist:
            session.add(TB.DRUG2SECEFF(words[3].lower(),words[6]))
        
        #And now just load them into the sqlite database 

## Ok, before we load it, we need to clean it and connect it with the ids we see
#  in the secondary effect database:

#Human_proteome_Loader()
#SecEff_Drugs_Scan()
#
## TODO: incorporate the percentages of the occurence

#for name in finlist:
#    session.add(TB.DrugName(name))

#SecEff_Drugs_Scan()
#SecEffIDs=list(set(SecEffIDs))
#
#for name in SecEffIDs:
#    session.add(TB.MedDRASecEffTerms(name[0],name[1]))

#session.commit()