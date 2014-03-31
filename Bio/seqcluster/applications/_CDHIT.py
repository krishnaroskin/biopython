# Copyright 2014 by Krishna M. Roskin.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the sequence clustering program CD-HIT.
"""

from __future__ import print_function

__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline

class CDHITCommandline(AbstractCommandline):
    """
    """
    def __init__(self, cmd="cd-hit", **kwargs):
        self.parameters = \
           [
            _Option(["-i", "i", "inputfile"],
                    "input filename in fasta format",
                    filename=True,
                    equate=False,
                    is_required=True),
            _Option(["-o", "o", "outputfile"],
                    "output filename",
                    filename=True,
                    equate=False,
                    is_required=True),
            _Option(["-c", "c", "cutoff"],
                    "sequence identity threshold, default 0.9",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, float) or isinstance(x, int)) and x >= 0.0 and x <= 1.0,
                    equate=False),
            _Option(["-G", "G", "global_ident"],
                    "use global sequence identity, default 1, if set to 0, then use local sequence identity",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and (x == 0 or x == 1),
                    equate=False),
            _Option(["-b", "b", "band_width"],
                    "band_width of alignment",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-M", "M", "memory"],
                    "memory limit (in MB) for the program, default 800; 0 for unlimitted",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-T", "T", "threads"],
                    "number of threads, default 1; with 0, all CPUs will be used",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-n", "n", "word_length"],
                    "word length",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-l", "l", "throw_len"],
                    "length of throw_away_sequences",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-t", "t", "redun_tol"],
                    "tolerance for redundance",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-d", "d", "des_len"],
                    "length of description in .clstr file, if set to 0, it takes the fasta defline and stops at first space",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-s", "s", "len_diff_cutoff"],
                    "length difference cutoff",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, int) or isinstance(x, float)) and x >= 0.0 and x <= 1.0,
                    equate=False),
            _Option(["-S", "S", "len_diff_cutoff_aa"],
                    "length difference cutoff in amino acid",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-aL", "aL", "align_cover_long"],
                    "alignment coverage for the longer sequence",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, int) or isinstance(x, float)) and x >= 0.0 and x <= 1.0,
                    equate=False),
            _Option(["-AL", "AL", "align_cover_long_aa"],
                    "alignment coverage control for the longer sequence",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-aS", "aS", "align_cover_short"],
                    "alignment coverage for the shorter sequence",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, int) or isinstance(x, float)) and x >= 0.0 and x <= 1.0,
                    equate=False),
            _Option(["-AS", "AS", "align_cover_short_aa"],
                    "alignment coverage control for the shorter sequence",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-A", "A", "align_cover"],
                    "minimal alignment coverage control for the both sequences",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-uL", "uL", "max_umatched_longer"],
                    "maximum unmatched percentage for the longer sequence",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, int) or isinstance(x, float)) and x >= 0.0 and x <= 1.0,
                    equate=False),
            _Option(["-uS", "uS", "max_umatched_shorter"],
                    "maximum unmatched percentage for the shorter sequence",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, int) or isinstance(x, float)) and x >= 0.0 and x <= 1.0,
                    equate=False),
            _Option(["-U", "U", "max_unmatch_len"],
                    "maximum unmatched length",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
            _Option(["-B", "B", "ram"],
                    "1 or 0, default 0, by default, sequences are stored in RAM, if set to 1, sequence are stored on hard drive",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and (x == 0 or x == 1),
                    equate=False),
            _Option(["-p", "p", "print_overlap"],
                    "default 0, if set to 1, print alignment overlap in .clstr file",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and (x == 0 or x == 1),
                    equate=False),
            _Option(["-g", "g", "global_search"],
                    "1 or 0, default 0, if set to 1, the program will cluster it into the most similar cluster that meet the threshold",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and (x == 0 or x == 1),
                    equate=False),
            ]
        AbstractCommandline.__init__(self, cmd, **kwargs)

def _test():
    """Run the module's doctests (PRIVATE)."""
    print("Running modules doctests...")
    import doctest
    doctest.testmod()
    print("Done")

if __name__ == "__main__":
    _test()
