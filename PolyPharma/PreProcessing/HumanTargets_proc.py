"""

"""
__author__ = 'ank'


infile='/home/akucahravy/Downloads/interaction.dat'
outfile='/home/akucahravy/Downloads/interaction_corr.dat'

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
