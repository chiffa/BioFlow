"""
Generic Setup configured for the project
"""
from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path
import os
from bioflow import __version__, __author__, __author_mail__

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'legacy_readme.md'), encoding='utf-8') as f:
    long_description = f.read()

on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

if not on_rtd:
    requirements = [
        'numpy',
        'scipy',
        'matplotlib',
        'scikit-learn',
        'cython',
        'pymongo',
        'requests',
        'click',
        'scikits.sparse',
        'mock',
        'requests-ftp',
        'neo4j-driver',
        'tabulate',
        'pyyaml',
    ]

else:
    requirements = [
        'pymongo',
        'requests',
        'click',
        'requests-ftp',
        'neo4j-driver',
    ]

setup(
    name='BioFlow',
    version=__version__,
    description='Information Flow Analysis in Biological Networks',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/chiffa/BioFlow',
    author=__author__,
    author_email=__author_mail__,
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
        'Programming Language :: Python :: 3.8',
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
               'coveralls<5.0',
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
