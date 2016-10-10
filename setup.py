#!/usr/bin/env python

from setuptools import setup
import os

VERSIONFILE=os.path.join("sixgill","version.py")
exec(open(VERSIONFILE).read())

long_description = 'six-frame genome-inferred libraries for LC-MS/MS'
if os.path.exists('README.rst'):
    long_description = open('README.rst').read()

version = open(VERSIONFILE, "rt").read().strip()

setup(name='sixgill',
      version=__version__,
      description='six-frame genome-inferred libraries for LC-MS/MS',
      author='Damon May',
      author_email='damonmay@uw.edu',
      packages=['sixgill'],
      license='Apache',
      install_requires=['biopython','pysam'],
      scripts=['bin/sixgill_build',
               'bin/sixgill_filter',
               'bin/sixgill_makefasta',
               'bin/sixgill_merge'],
      long_description=long_description,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: Apache Software License',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          ],
      keywords='proteomics metaproteomics metagenomics',
     )
