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
        'numpy',
        'scipy',
        'matplotlib',
        'scikit-learn',
        'python-Levenshtein',
        'cython',
        'bulbs',
        'pymongo',
        'requests',
        'click',
        'scikits.sparse',
        'mock',
        'requests-ftp']

else:
    requirements = [
        'bulbs',
        'pymongo',
        'click',
        'requests',
        'click',
        'mock',
        'requests-ftp']

setup(
    name='BioFlow',
    version='0.0.4',
    description='Information Flow Analyzer for biological networks',
    long_description=long_description,
    url='https://github.com/chiffa/BioFlow',
    author='Andrei Kucharavy',
    author_email='andrei.chiffa136@gmail.com',
    license='BSD',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research, Developers',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information analysis',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Operating System :: POSIX :: Linux',
        'License :: OSI Approved :: BSD 3-clause license',
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
