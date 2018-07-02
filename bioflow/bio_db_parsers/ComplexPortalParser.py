from csv import reader as csv_reader


def parse_complex_portal(complex_portal_file):

    def unpack_complex_contents(complex_name):
        unpacked_subnodes = []
        subnode_list = new_nodes[complex_name]
        for sub_node in subnode_list:
            if sub_node in complexes:
                unpacked_subnodes += unpack_complex_contents(sub_node)
            else:
                unpacked_subnodes +=sub_node

        return unpacked_subnodes


    complexes = []
    base = []
    new_nodes = []

    with open(complex_portal_file, 'rb') as source:
        reader = csv_reader(source, delimiter='\t')
        header = reader.next()
        for line in reader:
            legacy_id = line[0]
            display_name = line[1]
            componenets = line[4].split('|')
            componenets = [comp[:-3] for comp in componenets]
            node = {'ID': legacy_id, 'displayName': display_name, 'components': componenets}
            new_nodes.append(node)
            complexes.append(node['ID'])

    for node in new_nodes:
        node['components'] = unpack_complex_contents(node['ID'])
        base += node['components']

    base = list(set(base))

    return base, new_nodes

if __name__ == "__main__":
    parse_complex_portal('/home/andrei/sources/ComplexPortal/homo_sapiens.tsv')