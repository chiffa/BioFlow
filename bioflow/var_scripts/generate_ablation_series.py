from pathlib import Path
from csv import reader as csv_reader
from csv import writer as csv_writer
import random
import numpy as np

source_file = "C:\\Users\\Andrei\\Dropbox\\workspaces\\JHU\\Ewald Lab\\" \
              "Kp_Km data\\mouse_weighted_abs_log-fold.txt"

background_file = "C:\\Users\\Andrei\\Dropbox\\workspaces\\JHU\\Ewald Lab\\" \
                  "Kp_Km data\\mouse_genes_background.txt"

source_path = Path(source_file)

fname = source_path.stem
fname = fname + '_ablations'
storage_folder = Path(source_path.parent.joinpath(fname))
storage_folder.mkdir(parents=True, exist_ok=True)


def dump_experiment(experiment_name, corrected_lines):
    updated_experiment_name = experiment_name + '.tsv'
    print('debug', storage_folder)
    print('debug2', updated_experiment_name)
    with open(Path(storage_folder).joinpath(updated_experiment_name), 'wt',
              newline='', encoding='utf-8') as destination:
        writer = csv_writer(destination, delimiter='\t')
        writer.writerows(corrected_lines)


lines = []
with open(source_file, 'rt') as source:
    reader = csv_reader(source, delimiter='\t')
    for line in reader:
        lines.append(line)


background_lines = []
with open(background_file, 'rt') as source:
    reader = csv_reader(source, delimiter='\t')
    for line in reader:
        background_lines.append(line[0])


for removal_value in [0.05, 0.1, 0.2, 0.5]:

    # padding filter vector
    filter_line = [True]*len(lines)
    start_point = int(len(filter_line)*removal_value)
    filter_line = np.array(filter_line)
    filter_line[-start_point:] = False
    filter_line = filter_line.tolist()

    corrected_lines = [duplet for _filtered, duplet in zip(filter_line, lines) if _filtered]

    dump_experiment('lowest_%d_percent_removed' % (removal_value*100), corrected_lines)

    corrected_lines = [duplet
                       if _filtered
                       else [background_lines[random.randint(0, len(background_lines))], duplet[1]]
                       for _filtered, duplet
                       in zip(filter_line, lines)]

    dump_experiment('lowest_%d_percent_set_to_random' % (removal_value*100), corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [duplet for _filtered, duplet in zip(filter_line, lines) if _filtered]

    dump_experiment('random_%d_percent_removed' % (removal_value*100), corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [duplet
                       if _filtered
                       else [background_lines[random.randint(0, len(background_lines))], duplet[1]]
                       for _filtered, duplet
                       in zip(filter_line, lines)]

    dump_experiment('random_%d_percent_set_to_random' % (removal_value*100), corrected_lines)


flat_line = [(duplet[0],) for duplet in lines]
dump_experiment('no_weights', flat_line)

for removal_value in [0.05, 0.1, 0.2, 0.5]:

    # padding filter vector
    filter_line = [True]*len(lines)
    start_point = int(len(filter_line)*removal_value)
    filter_line = np.array(filter_line)
    filter_line[-start_point:] = False
    filter_line = filter_line.tolist()

    corrected_lines = [(duplet[0],) for _filtered, duplet in zip(filter_line, lines) if _filtered]
    dump_experiment('no_weights_lowest_%d_percent_removed' % (removal_value*100), corrected_lines)

    corrected_lines = [(duplet[0],)
                       if _filtered
                       else (background_lines[random.randint(0, len(background_lines))], )
                       for _filtered, duplet
                       in zip(filter_line, lines)]

    dump_experiment('no_weights_lowest_%d_percent_set_to_random' % (removal_value*100), corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [(duplet[0],) for _filtered, duplet in zip(filter_line, lines) if _filtered]
    dump_experiment('no_weights_random_%d_percent_removed' % (removal_value*100), corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [(duplet[0],)
                       if _filtered
                       else (background_lines[random.randint(0, len(background_lines))], )
                       for _filtered, duplet
                       in zip(filter_line, lines)]

    dump_experiment('no_weights_random_%d_percent_set_to_random' % (removal_value*100), corrected_lines)

# assumes that everything is in the order

