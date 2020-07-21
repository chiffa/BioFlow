from csv import reader as csv_reader


def parse_complex_portal(complex_portal_file):

    def unpack_complex_contents(complex_name):
        unpacked_subnodes = []
        subnode_list = new_nodes[complex_name]['components']

        for sub_node in subnode_list:

            if sub_node in list(new_nodes[complex_name].keys()):
                unpacked_subnodes += unpack_complex_contents(sub_node)

            else:
                if ':' in sub_node or '_9606' in sub_node:
                    pass
                elif '-' in sub_node:
                    unpacked_subnodes.append(sub_node.split('-')[0])
                else:
                    unpacked_subnodes.append(sub_node)

        return unpacked_subnodes

    base = []
    new_nodes = {}

    with open(complex_portal_file, 'rt') as source:
        reader = csv_reader(source, delimiter='\t')
        header = next(reader)
        for line in reader:
            legacy_id = line[0]
            display_name = line[1]
            componenets = line[4].split('|')
            componenets = [comp.split('(')[0] for comp in componenets]
            node = {'ID': legacy_id, 'displayName': display_name, 'components': componenets}
            new_nodes[node['ID']] = node

    # print new_nodes

    for node in new_nodes.values():
        node['components'] = unpack_complex_contents(node['ID'])
        base += node['components']

    # print new_nodes

    base = list(set(base))

    return new_nodes, base