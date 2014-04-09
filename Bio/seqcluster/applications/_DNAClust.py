# Copyright 2014 by Krishna M. Roskin.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Command line wrapper for the sequence clustering program DNACLUST.
"""

from __future__ import print_function

__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline

class DNAClustCommandline(AbstractCommandline):
    """Command line wrapper for the equence clustering program DNACLUST.

    http://dnaclust.sourceforge.net/

    Example:

    To cluster a FASTA file (seq.fasta) with DNACLUST, use:

    >>> from Bio.seqcluster.applications import DNAClustCommandline
    >>> dnaclust_cline = DNAClustCommandline(inputfile="seq.fasta",
    ...                                      similarity=0.45, header=True, threads=4)
    >>> print(dnaclust_cline)
    dnaclust --similarity 0.45 --header --threads 4 seq.fasta

    You would typically run the command line with dnaclust_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.
    Note that DNACLUST will write the alignment to stdout, which you may
    want to save to a file and then parse, e.g.::
    
        stdout, stderr = dnaclust_cline()
        with open("aligned.aln", "w") as handle:
            handle.write(stdout)

    Citations:

    Mohammadreza Ghodsi, Bo Liu and Mihai Pop. 2011. DNACLUST: accurate
    and efficient clustering of phylogenetic marker genes.
    BMC Bioinformatics 12: 271 

    Last checked against revision 54 dated: 2012-10-21 18:06:23 -0700 (Sun, 21 Oct 2012)
    """
    def __init__(self, cmd="dnaclust", **kwargs):
        self.parameters = \
           [
            _Option(["--similarity", "s", "similarity"],
                    "similarity threshold, between 0 and 1, specifing cluster radius",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, float) or isinstance(x, int)) and x >= 0.0 and x <= 1.0,
                    equate=False),
            _Option(["--predetermined-cluster-centers", "p", "predetermined_cluster_centers"],
                    "file containing predetermined cluster centers",
                    filename=True,
                    equate=False),
            _Switch(["--recruit-only", "r", "recruit_only"],
                    "cluster the input sequences to the predetermined centers "
                    "(used with --predetermined-cluster-centers)"),
            _Switch(["--header", "d", "header"],
                    "output header line with run options"),
            _Switch(["--left-gaps-allowed", "l", "left_gaps_allowed"],
                    "allow gaps on the left of shorter string in alignment"),
            _Option(["--k-mer-length", "k", "k_mer_length"],
                    "Maximum length of the k-mers used for filtering.",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
            _Switch(["--approximate-filter", "approximate_filter"],
                    "Use fast but approximate k-mer filter"),
            _Switch(["--no-k-mer-filter", "no_k_mer_filter"],
                    "disable k-mer filter"),
            _Switch(["--no-overlap", "no_overlap"],
                    "clusters are not allowed to overlap"),
            _Option(["--threads", "t", "threads"],
                    "number of threads to use",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
            _Switch(["--use-full-query-header", "u", "usefullqueryheader"],
                    "used the full sequence header instead of the first word"),
            _Argument(["input_file"],
                    "input file in FASTA",
                    filename=True,
                    is_required=True),
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
