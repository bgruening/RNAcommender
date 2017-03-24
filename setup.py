#!/usr/bin/env python
"""Setup RNAcommender."""

from setuptools import setup

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"
__version__ = 0.1


setup(
    name='rnacommender',
    version=__version__,
    author='Gianluca Corrado',
    author_email='gianluca.corrado@unitn.it',
    packages=['rnacommender',
              'rnacommender.utils',
              'rnacommender.fasta_utils',
              'rnacommender.pfam_utils'
              ],
    scripts=[],
    include_package_data=True,
    package_data={},
    license="MIT",
    description="""Genome-wide recommendation of RNA-protein interactions.""",
    long_description=open('README.md').read(),
    install_requires=[],
)
