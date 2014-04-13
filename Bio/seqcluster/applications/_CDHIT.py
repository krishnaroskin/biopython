# Copyright 2014 by Krishna M. Roskin.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Commandline wrappers for the family of sequence clustering programs CD-HIT.
"""

from __future__ import print_function

__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

from Bio.Application import _Option, _Switch, _Argument, AbstractCommandline

class _CDHITCommandlineBase(AbstractCommandline):
    """Base class for the all CD-HIT tools.

    This class serves as a base for all CD-HIT tools. Options shared between all CD-HIT
    tools are defined here.
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
                _Option(["-i", "i", "input_file"],
                    "input filename in fasta format",
                    filename=True,
                    equate=False,
                    is_required=True),
                _Option(["-o", "o", "output_file"],
                    "output filename",
                    filename=True,
                    equate=False,
                    is_required=True),
                _Option(["-c", "c", "cutoff"],
                    "sequence identity threshold",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, float) or isinstance(x, int)) and x >= 0.0 and x <= 1.0,
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
                _Option(["-aL", "aL", "long_align_coverage"],
                    "alignment coverage for the longer sequence",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, int) or isinstance(x, float)) and x >= 0.0 and x <= 1.0,
                    equate=False),
                _Option(["-AL", "AL", "long_align_coverage_aa"],
                    "alignment coverage control for the longer sequence",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
                _Option(["-aS", "aS", "short_align_coverage"],
                    "alignment coverage for the shorter sequence",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, int) or isinstance(x, float)) and x >= 0.0 and x <= 1.0,
                    equate=False),
                _Option(["-AS", "AS", "short_align_coverage_aa"],
                    "alignment coverage control for the shorter sequence",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
                _Option(["-B", "B", "ram"],
                    "1 or 0, default 0, by default, sequences are stored in RAM, if set to 1, sequence are stored on hard drive",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and (x == 0 or x == 1),
                    equate=False),
                _Option(["-g", "g", "global_search"],
                    "1 or 0, default 0, if set to 1, the program will cluster it into the most similar cluster that meet the threshold",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and (x == 0 or x == 1),
                    equate=False),
                ]
        # insert extra parameters at the start just in case there
        # are any arguments which must come last:
        self.parameters = extra_parameters + self.parameters
        AbstractCommandline.__init__(self, cmd, **kwargs)

class _CDHITOrCDHITESTCommandline(_CDHITCommandlineBase):
    """Base class for CD-HIT and CD-HIT-EST.

    This base class inherits from _CDHITCommandlineBase and add options shared with
    CD-HIT and CD-HIT-EST.
    """
    def __init__(self, cmd=None, **kwargs):
        assert cmd is not None
        extra_parameters = [
                _Option(["-G", "G", "global_identity"],
                    "use global sequence identity, default 1, if set to 0, then use local sequence identity",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and (x == 0 or x == 1),
                    equate=False),
                _Option(["-l", "l", "throw_away_len"],
                    "length of throw_away_sequences",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
                _Option(["-d", "d", "description_len"],
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
                _Option(["-A", "A", "min_align_coverage"],
                    "minimal alignment coverage control for the both sequences",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
                _Option(["-uL", "uL", "long_max_unmatched"],
                    "maximum unmatched percentage for the longer sequence",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, int) or isinstance(x, float)) and x >= 0.0 and x <= 1.0,
                    equate=False),
                _Option(["-uS", "uS", "short_max_unmatched"],
                    "maximum unmatched percentage for the shorter sequence",
                    filename=False,
                    checker_function=lambda x: (isinstance(x, int) or isinstance(x, float)) and x >= 0.0 and x <= 1.0,
                    equate=False),
                _Option(["-U", "U", "max_unmatch_len"],
                    "maximum unmatched length",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
                _Option(["-p", "p", "print_overlap"],
                    "default 0, if set to 1, print alignment overlap in .clstr file",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and (x == 0 or x == 1),
                    equate=False),
                ]
        # insert extra parameters at the start just in case there
        # are any arguments which must come last:
        self.parameters = extra_parameters + self.parameters
        _CDHITCommandlineBase.__init__(self, cmd, **kwargs)

class CDHITCommandline(_CDHITOrCDHITESTCommandline):
    """Command line wrapper for the sequence clustering program CD-HIT.

    http://weizhong-lab.ucsd.edu/cd-hit/

    Example:

    To cluster a FASTA file (seq.fasta) of amino acid sequences with CD-HIT, use:

    >>> from Bio.seqcluster.applications import CDHITCommandline
    >>> cdhit_cline = CDHITCommandline(input_file="seq.fasta", output_file="cluster",
    ...                                cutoff=0.80, print_overlap=1, threads=4)
    >>> print(cdhit_cline)
    cd-hit -i seq.fasta -o cluster -c 0.8 -T 4 -p 1

    You would typically run the command line with cdhit_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    The output file (cluster.clstr) can be parsed with Bio.seqcluster.CDHITClustIterator.

    Citations:

    "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide
    sequences", Weizhong Li & Adam Godzik Bioinformatics, (2006) 22:1658-9.

    Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu and Weizhong Li, CD-HIT: accelerated for
    clustering the next generation sequencing data. Bioinformatics, (2012),
    28 (23):3150-3152. doi: 10.1093/bioinformatics/bts565.

    Last checked against vision 4.5.4 dated: 2012-08-27.
    """
    def __init__(self, cmd="cd-hit", **kwargs):
        self.parameters = [
                _Option(["-t", "t", "redundance_tolerance"],
                    "tolerance for redundance",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and x >= 0,
                    equate=False),
                ]

        _CDHITOrCDHITESTCommandline.__init__(self, cmd, **kwargs)

class CDHITESTCommandline(_CDHITOrCDHITESTCommandline):
    """Command line wrapper for the sequence clustering program CD-HIT-EST.

    http://weizhong-lab.ucsd.edu/cd-hit/

    Example:

    To cluster a FASTA file (seq.fasta) of nucleotide sequences with CD-HIT-EST, use:

    >>> from Bio.seqcluster.applications import CDHITESTCommandline
    >>> cdhitest_cline = CDHITESTCommandline(input_file="seq.fasta", output_file="cluster",
    ...                                      cutoff=0.80, threads=4)
    >>> print(cdhitest_cline)
    cd-hit-est -i seq.fasta -o cluster -c 0.8 -T 4

    You would typically run the command line with cdhitest_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    The output file (cluster.clstr) can be parsed with Bio.seqcluster.CDHITClustIterator.

    Citations:

    "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide
    sequences", Weizhong Li & Adam Godzik Bioinformatics, (2006) 22:1658-9.

    Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu and Weizhong Li, CD-HIT: accelerated for
    clustering the next generation sequencing data. Bioinformatics, (2012),
    28 (23):3150-3152. doi: 10.1093/bioinformatics/bts565.

    Last checked against vision 4.5.4 dated: 2012-08-27.
    """
    def __init__(self, cmd="cd-hit-est", **kwargs):
        self.parameters = [
                _Option(["-r", "r", "search_reverse_comp"],
                    "search on both strands if 1, the default, if 0 only search on positive strand",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int) and (x == 0 or x == 1),
                    equate=False),
                _Option(["-match", "match", "matching_score"],
                    "matching score, default 2 (1 for T-U and N-N)",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
                _Option(["-mismatch", "mismatch", "mismatching_score"],
                    "mismatching score, default -2",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
                _Option(["-gap", "gap", "gap_open_score"],
                    "gap opening score, default -6",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
                _Option(["-gap-ext", "gap_ext", "gap_extend_score"],
                    "gap extension score, default -1",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
                ]

        _CDHITOrCDHITESTCommandline.__init__(self, cmd, **kwargs)

class CDHIT454Commandline(_CDHITCommandlineBase):
    """Command line wrapper for the sequence clustering program CD-HIT-454.

    http://weizhong-lab.ucsd.edu/cd-hit/

    Example:

    To cluster a FASTA file (seq.fasta) of nucleotide sequences with CD-HIT-454, use:

    >>> from Bio.seqcluster.applications import CDHIT454Commandline
    >>> cdhit454_cline = CDHIT454Commandline(input_file="seq.fasta", output_file="cluster",
    ...                                      cutoff=0.80, threads=4)
    >>> print(cdhit454_cline)
    cd-hit-454 -i seq.fasta -o cluster -c 0.8 -T 4

    You would typically run the command line with cdhit454_cline() or via
    the Python subprocess module, as described in the Biopython tutorial.

    The output file (cluster.clstr) can be parsed with Bio.seqcluster.CDHITClustIterator.

    Citations:

    "Cd-hit: a fast program for clustering and comparing large sets of protein or nucleotide
    sequences", Weizhong Li & Adam Godzik Bioinformatics, (2006) 22:1658-9.

    Limin Fu, Beifang Niu, Zhengwei Zhu, Sitao Wu and Weizhong Li, CD-HIT: accelerated for
    clustering the next generation sequencing data. Bioinformatics, (2012),
    28 (23):3150-3152. doi: 10.1093/bioinformatics/bts565.

    Last checked against vision 4.5.4 dated: 2012-08-27.
    """
    def __init__(self, cmd="cd-hit-454", **kwargs):
        self.parameters = [
                _Option(["-D", "D", "max_indel_size"],
                    "max size per indel, default 1",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
                _Option(["-match", "match", "matching_score"],
                    "matching score, default 2",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
                _Option(["-mismatch", "mismatch", "mismatching_score"],
                    "mismatching score, default -1",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
                _Option(["-gap", "gap", "gap_open_score"],
                    "gap opening score, default -3",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
                _Option(["-gap-ext", "gap_ext", "gap_extend_score"],
                    "gap extension score, default -1",
                    filename=False,
                    checker_function=lambda x: isinstance(x, int),
                    equate=False),
                ]
        _CDHITCommandlineBase.__init__(self, cmd, **kwargs)

def _test():
    """Run the module's doctests (PRIVATE)."""
    print("Running modules doctests...")
    import doctest
    doctest.testmod()
    print("Done")

if __name__ == "__main__":
    _test()
