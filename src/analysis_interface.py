"""
Provides the implementation of routines responsible for the analysis implementation
"""
from src.annotation_network.BioKnowledgeInterface import GeneOntologyInterface as AnnotomeInterface, get_background
from src.annotation_network.knowledge_access_analysis import auto_analyze as knowledge_analysis
from src.molecular_network.InteractomeInterface import InteractomeInterface as InteractomeInterface
from src.molecular_network.interactome_analysis import auto_analyze as interactome_analysis


def analyze_interactome(source, background=None, depth=100, processors=2):
    local_matrix = InteractomeInterface(main_connex_only=True, full_impact=True)
    local_matrix.full_rebuild()
    if background is not None:
        background = get_background(background)
    source_set = get_background(source)
    interactome_analysis(source_set, depth, processors, background)
    print "analsysis is finished, current results are stored in the $PROJECT_HOME/src/outputs directory"


def analyze_annotome(source, background=None, depth=100, processors=2):
    local_matrix = InteractomeInterface(main_connex_only=True, full_impact=True)
    local_matrix.full_rebuild()
    filtr = ['biological_process']
    if background is None:
        annot_matrix = AnnotomeInterface(
            filtr, local_matrix.all_uniprots_bulbs_id_list, (1, 1), True, 3)
    else:
        background_set = get_background(background)
        annot_matrix = AnnotomeInterface(
            filtr, background_set, (1, 1), True, 3)
    annot_matrix.rebuild()
    annot_matrix.store()
    source_set = get_background(source)
    knowledge_analysis(
        source=source_set,
        KG_object=annot_matrix,
        desired_depth=depth,
        processors=processors)
    print "analsysis is finished, current results are stored in the $PROJECT_HOME/src/outputs directory"

if __name__ == "__main__":
    # analyze_interactome()
    # analyze_annotome()
    pass