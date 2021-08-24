from pathlib import Path
from csv import reader as csv_reader
from csv import writer as csv_writer
import random

source_file = ''
background_file = ''

source_path = Path(source_file)

fname = source_path.stem
fname = fname + '_ablations'
storage_folder = Path(source_path.parent.joinpath(fname)).mkdir(parents=True, exist_ok=True)


def dump_experiment(experiment_name, corrected_lines):
    with open(Path(storage_folder).joinpath(experiment_name + '.tsv'), 'wt') as destination:
        writer = csv_writer(destination, delimiter='\t')
        writer.writerows(corrected_lines)


lines = []
with open(source_file, 'rt') as source:
    reader = csv_reader(source, delimiter='\t')
    for line in reader:
        lines.append(lines)


background_lines = []
with open(background_file, 'rt') as source:
    reader = csv_reader(source, delimiter='\t')
    for line in reader:
        background_lines.append(line[0]) # INTEST: make sure it is parsed properly

    print(background_lines)


for removal_value in [0.05, 0.1, 0.2, 0.5]:
    filter_line = [True]*len(lines)
    filter_line[int(len(filter_line)*removal_value)] = False

    corrected_lines = [duplet for _filtered, duplet in zip(filter_line, lines) if _filtered]
    dump_experiment('lowest_%d_percent_removed' % removal_value*100, corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [duplet
                       if _filtered
                       else [background_lines[random.randint(0, len(background_lines))], duplet[0]]
                       for _filtered, duplet
                       in zip(filter_line, lines) ]

    dump_experiment('lowest_%d_percent_set_to_random' % removal_value*100, corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [duplet for _filtered, duplet in zip(filter_line, lines) if _filtered]
    dump_experiment('random_%d_percent_removed' % removal_value*100, corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [duplet
                       if _filtered
                       else [background_lines[random.randint(0, len(background_lines))], duplet[0]]
                       for _filtered, duplet
                       in zip(filter_line, lines) ]

    dump_experiment('random_%d_percent_set_to_random' % removal_value*100, corrected_lines)


flat_line = [duplet[0] for duplet in lines]
dump_experiment('no_weights', flat_line)

for removal_value in [0.05, 0.1, 0.2, 0.5]:
    filter_line = [True]*len(lines)
    filter_line[int(len(filter_line)*removal_value)] = False

    corrected_lines = [duplet[0] for _filtered, duplet in zip(filter_line, lines) if _filtered]
    dump_experiment('no_weights_lowest_%d_percent_removed' % removal_value*100, corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [duplet[0]
                       if _filtered
                       else background_lines[random.randint(0, len(background_lines))]
                       for _filtered, duplet
                       in zip(filter_line, lines) ]

    dump_experiment('no_weights_lowest_%d_percent_set_to_random' % removal_value*100, corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [duplet[0] for _filtered, duplet in zip(filter_line, lines) if _filtered]
    dump_experiment('no_weights_random_%d_percent_removed' % removal_value*100, corrected_lines)

    random.shuffle(filter_line)
    corrected_lines = [duplet[0]
                       if _filtered
                       else background_lines[random.randint(0, len(background_lines))]
                       for _filtered, duplet
                       in zip(filter_line, lines) ]

    dump_experiment('no_weights_random_%d_percent_set_to_random' % removal_value*100, corrected_lines)

# assumes that everything is in the order

