#!/usr/bin/env python
"""
Register sixgill with PyPI, first converting README.md to .rst
"""

#import pandoc
import os
import sys

#pandoc.core.PANDOC_PATH = '/usr/local/bin/pandoc'

#doc = pandoc.Document()
#doc.markdown = open('README.md').read()
#f = open('README.txt', 'w+')
#f.write(doc.rst)
#f.close()
command = "python setup.py register"
if len(sys.argv) > 1:
    if sys.argv[1] == 'test':
        command = command + " -r https://testpypi.python.org/pypi"
    else:
        quit("Bad args.")
print(command)
os.system(command)
#os.remove('README.txt')
