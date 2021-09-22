"""
Records the run metadata.
The underlying object is a .yaml file containg a stack of all the records

[Downloads]
    - DB1:
        - download date
        - hash
        - download location
    - ...
[neo4jbuild]
    - build date
    - organism (check consistency with the current one)
[laplaciansbuild]
    - build date
[translation to internal IDs]
    - primary file that is being translated
    - secondary file that is being translated
    - tertiary file that is being translatede
[code]
    - current code version
    - python version
    - dependent libraries versions
    - current git commit shorthash
    - if any modifications are present  outside the analysis_pipeline_example [WARN the user to
    create a commit before running any real analysis]

inside the .meta folder in the run file, copy
    - Create  a copy of the active .confs file
    - Create a copy of the active data record stack
    - Create a copy of the
    - create a registration for every translation step.
"""
import yaml
from bioflow import __version__


class RecordsStack(object):

    def __init__(self, expected_location):
        """

        :param expected_location: where the record is supposed to be located on the drive
        """
        # create if doesn't exist:
        pass

        # read:
        with open(expected_location, 'r', encoding='utf8') as records_stack_f:
            records_stack = yaml.safe_load(records_stack_f)

        self.records_stack = records_stack

    def organisms_are_compatible(self, organism):
        """
        checks if the organism against wthich the db and laplacians were build are the same as
        currently active

        :param organism: current organism
        :return:
        """
        try:
            org = self.records_stack['neo4jbuild']['organism']
        except KeyError:
            raise Exception('organism undefined')

        if org != organism:
            raise Exception('different organisms in configs and ')

    def add_to_run(self):
        """
        Adds the current

        :return:
        """
        pass