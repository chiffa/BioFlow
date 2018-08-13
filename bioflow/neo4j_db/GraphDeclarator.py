"""
An interface declaration required by bulbs module
"""
from bioflow.main_configs import neo4j_server
from bioflow.utils.log_behavior import get_logger
from cypher_drivers import GraphDBPipe
import os


log = get_logger(__name__)

on_rtd = os.environ.get('READTHEDOCS') == 'True'
on_unittest = os.environ.get('UNITTESTING') == 'True'


if on_rtd or on_unittest:
    log.debug(
        'graph database interface mocks DB connection instead of actually connecting to it')

    from mock import Mock as MagicMock

    class Mock(MagicMock):

        @classmethod
        def __getattr__(cls, name):
            return Mock()

        @classmethod
        def __getitem__(cls, name):
            return Mock()

    DatabaseGraph = Mock()

else:
    DatabaseGraph = GraphDBPipe()
    DatabaseGraph.build_indexes()

