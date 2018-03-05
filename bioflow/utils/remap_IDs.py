"""
Remaps the IDs of the gene identifiers to a different organism before processing them
"""
from csv import reader as csv_reader
from csv import writer as csv_writer

high_conf_translation_dict = {}
low_conf_translation_dict = {}
genes_to_ids_dict = {}

# translation_file_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Veena data/Mouse_2_human.tsv'
# gene_to_id_file_location = ''
# data_source_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Veena data/both_ENSMUG.csv'
# data_dump_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Veena data/both_ENSHUM.csv'


translation_file_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Veena data/Mouse_2_human.tsv'
gene_to_id_file_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Kp_Km data/mouse_look_up_table.tsv'
data_source_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Kp_Km data/all_significant.csv'
data_dump_location = '/home/andrei/Dropbox/workspaces/JHU/Ewald Lab/Kp_Km data/all_sig_hum.csv'


with open(translation_file_location, 'r') as source:
    reader = csv_reader(source, delimiter='\t')
    print reader.next()
    for line in reader:
        if line[0] and line[1]:
            if int(line[3]):
                # We still need to account for the confidence in mapping
                high_conf_translation_dict[line[0]] = [line[1], line[2]]
                # print line[0:4]
            else:
                low_conf_translation_dict[line[0]] = [line[1], line[2]]

high_conf_trans = []
low_conf_trans = []


if gene_to_id_file_location:
    with open(gene_to_id_file_location, 'r') as source:
        reader = csv_reader(source, delimiter='\t')
        print reader.next()
        for line in reader:
            genes_to_ids_dict[line[2]] = line[0]


with open(data_source_location, 'r') as source:
    reader = csv_reader(source)
    for i, line in enumerate(reader):
        word = line[0]
        if gene_to_id_file_location:
            word = genes_to_ids_dict.get(word, 'None found')
        if word in high_conf_translation_dict.keys():
            high_conf_trans.append(high_conf_translation_dict[word])
        if word in low_conf_translation_dict.keys():
            low_conf_trans.append(low_conf_translation_dict[word])

print "out of %s, %s were translated with high confidence, %s with low and %s were not found" % \
      (i, len(high_conf_trans), len(low_conf_trans), i-len(high_conf_trans)-len(low_conf_trans))

with open(data_dump_location, 'w') as destination:
    writer = csv_writer(destination)
    writer.writerows((word for word in high_conf_trans))



