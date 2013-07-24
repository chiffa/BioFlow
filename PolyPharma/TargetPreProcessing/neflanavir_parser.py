'''
Created on Jul 23, 2013

@author: andrei
'''

subdict={
'cell division protein kinase 6' : 'cell division protein kinase 6',
'cell division protein kinase 2' : 'cell division protein kinase 2',
'reelin' : 'reelin',
'serine/threonine-protein kinase B-RAF' : 'serine/threonine-protein kinase b-raf',
'gastricsin' : 'gastricsin',
'serine/threonine protein kinase TAO2' : 'serine/threonine-protein kinase tao2',
'ephrin type-B receptor 4' : 'ephrin type-b receptor 4',
'beta-actin' : 'beta-actin',
'activated cdc42 kinase 1' : 'activated cdc42 kinase 1',
'cyclin-dependent kinase 4' : 'cyclin-dependent kinase 4',
'3-phosphoinositide dependent protein kinase 1' : '3-phosphoinositide-dependent protein kinase 1',
'pantothenate kinase 1' : 'pantothenate kinase 1',
'beta-secretase 1' : 'beta-secretase 1',
'beta-secretase 2' : 'beta-secretase 2',
'protein-l-isoaspartate(D-aspartate)-o-methyltransferase' : 'protein-l-isoaspartate(d-aspartate) o-methyltransferase',
'cyclin-dependent protein kinase 2' : 'cyclin-dependent kinase 2',
'glycerol kinase' : 'glycerol kinase',
'cyclin-dependent protein kinase 5' : 'cyclin-dependent kinase 5',
'protein kinase C, theta type' : 'protein kinase c theta type',
'tyrosine-protein kinase ABL1' : 'tyrosine-protein kinase abl1',
'tyrosine-protein kinase HCK' : 'tyrosine-protein kinase hck',
'mitogen-activated protein kinase 8' : 'mitogen-activated protein kinase 8',
'tyrosine-protein kinase lck' : 'tyrosine-protein kinase lck',
'epidermal growth factor receptor' : 'epidermal growth factor receptor',
'Rho-associated protein kinase 2' : 'rho-associated protein kinase 2',
'angiopoietin-1 receptor' : 'angiopoietin-1 receptor',
'serine/threonine kinase 6' : 'serine/threonine-protein kinase 6',
'm-calpain' : 'm-calpain',
'transitional endoplasmic reticulum atpase' : 'transitional endoplasmic reticulum atpase',
'serine/threonine-protein kinase PIM-1' : 'serine/threonine-protein kinase pim-1',
'G-protein coupled receptor kinase 2' : 'g-protein coupled receptor kinase 2',
'ephrin type-A receptor 2' : 'ephrin type-a receptor 2',
'insulin receptor subunit beta' : 'insulin receptor subunit beta',
'cytochrome b' : 'cytochrome b',
'interleukin-1 receptor-associated kinase 4' : 'interleukin-1 receptor-associated kinase 4',
'ephrin type-A receptor 7' : 'ephrin type-a receptor 7',
'ephrin type-A receptor 5' : 'ephrin type-a receptor 5',
'cytochrome p450 2c8' : 'cytochrome p450 2c8',
'renin' : 'renin',
'rho-associated protein kinase 1' : 'rho-associated protein kinase 1',
'tyrosine-protein kinase mer' : 'tyrosine-protein kinase mer',
'aurora-related kinase 2' : 'aurora-related kinase 2',
'focal adhesion kinase 1' : 'focal adhesion kinase 1',
'ADP-dependent glucokinase' : 'adp-dependent glucokinase',
'serine/threonine-protein kinase wnk1' : 'serine/threonine-protein kinase wnk1',
'vitamin d binding protein' : 'vitamin d-binding protein',
'prohibitin' : 'prohibitin',
'vascular endothelial growth factor receptor 2' : 'vascular endothelial growth factor receptor 2',
'fibroblast growth factor receptor 2' : 'fibroblast growth factor receptor 2', #Insure from here
'pyrroline-5-carboxylate reductase 1':'pyrroline-5-carboxylate reductase 1, mitochondrial',
'dual specificity mitogen-activated protein  kinase 1':'dual specificity mitogen-activated protein  kinase 1',
'RAC-beta serine/threonine-protein kinase (AKT2)':'rac-beta serine/threonine-protein kinase',
'cyclin-dependent protein kinase pho85':'', #Not a human protein
'tyrosine-protein kinase receptor RET':['RET_HUMAN'],
'adenylyltransferase THIF':'', #Not a human protein
'heat shock cognate':'heat shock cognate 71 kda protein',
'serine/threonine-protein kinase 24':['STK24_HUMAN'],
'erine/threonine protein kinase TAO2':'serine/threonine-protein kinase tao2',
'acetyl-coenzyme a synthetase':'acetyl-coa synthetase',
'tyrosine-protein kinase src':['SRC_HUMAN'],
'cAMP-dependent protein kinase':'', #??????
'amyloid protein-binding protein':'amyloid protein-binding protein 2',
'casein kinase 1 gamma 2':'casein kinase i isoform gamma-2',
'dihydropyrimidine dehydrogenase':'dihydropyrimidine dehydrogenase [nadp(+)]',
'DNA polymerase iii alpha subunit':'', #Not a human protein
'NADH-quinone oxidoreductase':[''], #Several genes (NAO*_Human)
'complement c3 beta chain':['CO3_HUMAN'],
'insulin-like growth factor 1 receptor':'insulin-like growth factor i receptor',
'serine/threonine-protein kinase pknB':'', #Not a human protein
'pyruvate dehydrogenase kinase': ['PDK1_HUMAN','PDK2_HUMAN','PDK3_HUMAN','PDK4_HUMAN'],
'kinesin heavy chain-like protein':'', #Not a human protein
'pepsin 3A':'pepsin a-3',
'nagk protein':['NAGK_HUMAN'],
'SAM-dependent methyltransferase':'', #Not a human protein
'tyrosine-protein kinase ERBB-4':['ERBB4_HUMAN'],
'pyruvate kinase':'pyruvate kinase 1',
'phosphatidylinositol 3-kinase':'', #Pretty huge group of kinases....
'TGF-A superfamily receptor type I':'', #Is it possible to get the refernces for the elements of this family?
'RNA uridylyl transferase':'', # Not a protein????
'putative asparaginyl hydroxylase':'aspartyl/asparaginyl beta-hydroxylase',
'acetylcholine receptor protein, gamma chain':'acetylcholine receptor subunit gamma',
'SAM-dependent methyltransferase, REBM':'', # Does not exist in Uniprot
'PSD-95 SH3-guanylate kinase domain':['DLG4_HUMAN'],    # 
'benzoate-coenzyme A ligase':'', # Not a human protein
'IGF-1 receptor kinase':['IGF1R_HUMAN'],
'ribosomal protein':'', #Totally unsure from here
'dihydroxyacetone kinase':['DHAK_HUMAN'],
'serine threonine kinase endoribonuclease ire1':['ERN1_HUMAN'],
'3 isopropylmalate dehydrogenase':'', # Does not seem to be human. Isocytrate?
'uropepsin':'', # Does not exists in Uniport. maybe it is urotensin?
'fgf receptor 1 kinase domain':['FGFR1_HUMAN'],
'ero1p':['ERO1A_HUMAN']
}
