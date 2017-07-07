#!/usr/bin/env sage
# -*- mode: sage -*-
r"""
Run doctests
"""

from sage.doctest.control import DocTestDefaults, DocTestController

# Add any new files to this list
files = [
    'hillman_grassl_tableau.py',
    'reverse_plane_partition.py'
]

DD = DocTestDefaults()
DC = DocTestController(DD, files)
DC.run()
