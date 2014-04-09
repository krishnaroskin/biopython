# Copyright 2014 by Krishna M. Roskin.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Code for dealing with sequence clustering.
"""
from __future__ import print_function

import StringIO
import tempfile
import shutil
import os.path

from Bio import SeqIO

from .applications import DNAClustCommandline
from .applications import CDHITESTCommandline

__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

class SeqCluster(object):
    def __init__(self, name, representative, members):
        self.name = name
        self.representative = representative
        self.members = members

    def __iter__(self):
        return self.members.__iter__()

class SeqClusterMember(object):
    def __init__(self, name, percent_match=None):
        self.name = name
        self.percent_match = percent_match

def DNAClustIterator(handle):
    """
    """
    # skip any header stuff
    line = handle.readline()
    while True:
        if line == "":
            return  # Premature end of file, or just empty?
        elif line.startswith("%"):  # skip lines comment lines
            line = handle.readline()
        else:
            break

    while True:
        # parse the cluster name, use the first/representative sequence name as the cluster name
        cluster_members = [SeqClusterMember(name) for name in line.strip().split("\t")]
        line = handle.readline()

        yield SeqCluster(cluster_members[0].name, cluster_members[0], cluster_members)

        if not line:    # end of file
            return

    assert False, "Should not reach this line"

def ParseCDHITMember(line):
    """
    """
    # strip off the cluster number
    number, rest = line.split("\t")
    number = int(number)

    # strip off the length
    length, rest = rest.split(", >", 1)
    assert length.endswith("aa") or length.endswith("nt"), "Unknown length units for %s" % length
    length = int(length[:-2])

    strand       = None
    rep_start    = None
    rep_stop     = None
    member_start = None
    member_stop  = None
    identity     = None

    if "... at " in rest:   # this is normal member
        # strip off the sequence name
        name, rest = rest.split("... at ", 1)

        # strip off the overlap and/or strand
        if rest.count("/") == 1:     # if ouput includes overlap or strand
            if rest.startswith("+") or rest.startswith("-"):    # figure out which
                strand, identity = rest.split("/")
            else:
                overlap, identity = rest.split("/")
                rep_start, rep_stop, member_start, member_stop = overlap.split(":")
                rep_start    = int(rep_start)
                rep_stop     = int(rep_stop)
                member_start = int(member_start)
                member_stop  = int(member_stop)
        elif rest.count("/") == 2:  # both strand and overlap are present
            overlap, strand, identity = rest.split("/")
            rep_start, rep_stop, member_start, member_stop = overlap.split(":")
            rep_start    = int(rep_start)
            rep_stop     = int(rep_stop)
            member_start = int(member_start)
            member_stop  = int(member_stop)
        else:
            identity = rest

        # parse the identity
        assert identity.endswith("%\n"), "Identity should end with '%' "
        identity = float(identity[:-len("%\n")])

    else:   # this is the representative member
        # strip off the sequence name
        assert rest.endswith("... *\n")
        name = rest[:-len("... *\n")]

    return number, length, name, (rep_start, rep_stop), (member_start, member_stop), identity, strand

def CDHITClustIterator(handle):
    """
    """
    # skip any header stuff, i.e. everything until we hit a >
    line = handle.readline()
    while True:
        if line == "":
            return  # premature end of file, or just empty?
        elif line.startswith(">"):  # new clusters start with >
            break

    while True:
        # new clusters start with >
        if line[0] != ">":
            raise ValueError("CD-HIT files should start with a '>' character")
        cluster_name = line[1:].strip()
        cluster_members = []

        line = handle.readline()
        while True:
            if not line:    # end of file
                break
            elif line[0] == ">":    # start of a new cluster
                break
            
            # parse the member line
            number, length, name, rep_overlap, member_overlap, identity, strand = ParseCDHITMember(line)
            new_member = SeqClusterMember(name, identity)
            # if this is the representative (indicated by percent identity being None), store it
            if identity is None:
                representative = new_member
            cluster_members.append(new_member)

            line = handle.readline()

        yield SeqCluster(cluster_name, representative, cluster_members)

        if not line:    # end of file
            return  # StopIteration

    assert False, "Should not reach this line"

def DNAClustHandler(sequences, identity_cutoff, **kwargs):

    # write the sequences to a temp. FASTA file
    tmp_file = tempfile.NamedTemporaryFile("w")
    SeqIO.write(sequences, tmp_file, "fasta")
    tmp_file.flush()

    cmd = DNAClustCommandline(input_file=tmp_file.name,
                              similarity=identity_cutoff,
                              **kwargs)
    
    stdout, stderr = cmd()
    return DNAClustIterator(StringIO.StringIO(stdout))

def CDHITEXTClustHandler(sequences, identity_cutoff, **kwargs):
    # make temp. directory
    temp_dir_path = tempfile.mkdtemp()
    import sys
    print(sys.stderr, temp_dir_path)

    try:
        sequence_filename = os.path.join(temp_dir_path, "seq.fasta")
        sequence_file = open(sequence_filename, "w")
        SeqIO.write(sequences, sequence_file, "fasta")
        sequence_file.close()

        cmd = CDHITESTCommandline(input_file=sequence_filename,
                                  output_file=os.path.join(temp_dir_path, "clusters"),
                                  description_len=0,
                                  **kwargs)
        cmd()
        
        # read all the clusters since we're going to delete everything
        clusters = list(CDHITClustIterator(open(os.path.join(temp_dir_path, "clusters.clstr"), "r")))
    finally:
        shutil.rmtree(temp_dir_path)

    return clusters

if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
