# Copyright 2014 by Krishna M. Roskin.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Alignment command line tool wrappers."""

__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

from ._DNAClust import DNAClustCommandline

#Make this explicit, then they show up in the API docs
__all__ = ["DNAClustCommandline"]
