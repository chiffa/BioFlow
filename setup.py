"""
Generic Setup configured for the project
"""
from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path
import os

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:
    requirements = [
        'numpy < 2.0',
        'scipy < 1.0',
        'matplotlib < 2.0',
        'scikit-learn < 0.17',
        'cython < 0.23',
        'pymongo < 4.0',
        'requests < 3.0',
        'click < 6.0',
        'scikits.sparse < 0.3',
        'mock < 2.0',
        'requests-ftp < 0.4',
        'neo4j-driver < 2.0',
    ]

else:
    requirements = [
        'pymongo < 4.0',
        'requests < 3.0',
        'click < 6.0',
        'requests-ftp < 0.4',
        'neo4j-driver < 2.0',]

setup(
    name='BioFlow',
    version='0.2.2',
    description='Information Flow Analysis in biological networks',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/chiffa/BioFlow',
    author='Andrei Kucharavy',
    author_email='andrei.chiffa136@gmail.com',
    license='BSD',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Operating System :: POSIX :: Linux',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    keywords='network analysis, systems biology, interactome, computational biology',
    packages=['bioflow']+find_packages(),
    package_data={'bioflow': ['configs/*.ini']},
    include_package_data=True,
    install_requires=requirements,
    extras_require={
        'dev': [],
        'test': [],
        'ci': ['coverage',
               'pylint',
               'coveralls',
               'pyflakes',
               'pep-8-naming',
               'maccabe',
               'flake-8'],
    },
    entry_points="""
    [console_scripts]
    bioflow = bioflow.cli:main
    """,
)
