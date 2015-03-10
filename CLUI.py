__author__ = 'ank'

import click


@click.group()
def truegird():
    pass




@click.command()
def verify():
    """
    Verifies that all has been properly installed
    :return:
    """
    pass


@click.command()
@click.argument('organism')  #help='name of the organism database to load. Current options: human, mouse, yeast')
def load_database(organism):
    pass


@click.command()
def fill_database():
    pass




truegird.add_command(verify)
truegird.add_command(fill_database)

if __name__ == '__main__':
    truegird()