#!/usr/bin/env sage
# -*- mode: sage -*-
r"""
Run doctests
"""

from sage.doctest.control import DocTestDefaults, DocTestController
import os
import sys

# Add any new files to this list
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

DD = DocTestDefaults()
DC = DocTestController(DD, ['ptdt_package', 'chern_char.sage'])
DC.run()
