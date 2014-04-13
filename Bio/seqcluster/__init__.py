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
from .applications import CDHIT454Commandline
from .applications import CDHITESTCommandline
from .applications import CDHITCommandline

__docformat__ = "epytext en"  # Don't just use plain text in epydoc API pages!

class SeqCluster(object):
    """A cluster of sequence based objects.

    Class to store a cluster of sequence objects, the cluster members, and the
    cluster representative.
    """
    def __init__(self, name, representative, members):
        self.name = name
        self.representative = representative
        self.members = members

    def __iter__(self):
        """
        Iterating over a cluster is to iterate over its members.
        """
        return self.members.__iter__()

class SeqClusterMember(object):
    """A member of a sequence based cluster.

    Only the name is stored in the cluster member. No sequene is store. That must be
    stored elsewhere.
    """
    def __init__(self, name):
        self.name = name

def DNAClustIterator(handle):
    """Parse the outout of dnaclust.

    Parsed the output of dnaclust and yeilds SeqCluster objects for each cluster. The name
    of the cluster is the name of the representative sequence, which is also the first
    sequence given.
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
        # parse the cluster members
        cluster_members = [SeqClusterMember(name) for name in line.strip().split("\t")]
        # fetch the next line
        line = handle.readline()

        # the cluster name is the name of the first/representative sequence
        yield SeqCluster(cluster_members[0].name, cluster_members[0], cluster_members)

        if not line:    # end of file
            return

    assert False, "Here there be dagons"

def ParseCDHITMember(line):
    """Parse CD-HIT member definition line.

    Parse and return information derived from a CD-HIT/CD-HIT-EST/CD-HIT-454 member
    definition.
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
        identity = float(identity[:-2])     # convert after removing % and \n

    else:   # this is the representative member
        # strip off the sequence name
        assert rest.endswith("... *\n")
        name = rest[:-6]    # remove the trailing ... *\n

    return number, length, name, (rep_start, rep_stop), (member_start, member_stop), identity, strand

def CDHITClustIterator(handle):
    """Parse the outout of one of the CD-HIT programs.

    Parsed the output of CD-HIT, CD-HIT-EST, or CD-HIT-454 and yeilds SeqCluster objects
    for each cluster. The name of the cluster is given as "Cluster NNN" as given by CD-HIT.
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
        assert line[0] == ">", "CD-HIT output files should start with a '>' character"
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
            new_member = SeqClusterMember(name)
            # add the CD-HIT specific metadata, identity handled below
            new_member.length = length
            if rep_overlap:
                assert member_overlap is not None
                new_member.rep_overlap = rep_overlap
                new_member.member_overlap = member_overlap
            if strand:
                new_member.strand = strand
            # if this is the representative (indicated by percent identity being None), store it
            if identity is None:
                representative = new_member
                new_member.identity = 100.0
            else:
                new_member.identity = identity

            cluster_members.append(new_member)

            line = handle.readline()    # get the next line

        yield SeqCluster(cluster_name, representative, cluster_members)

        if not line:    # end of file
            return  # StopIteration

    assert False, "I found a counter-example to the Riemann hypothesis!"

def DNAClustHandler(sequences, identity_cutoff, **kwargs):
    """Cluster DNA sequences with dnaclust.
    
    Clusters the given DNA sequences with dnaclust and return yield SeqCluster objects
    for each cluster. The cluster radius is given by identity_cutoff. Other keyword
    parameters are passed to dnaclust.
    """
    # TODO Check that input sequences are DNA

    # write the sequences to a temp. FASTA file
    tmp_file = tempfile.NamedTemporaryFile("w")
    SeqIO.write(sequences, tmp_file, "fasta")
    tmp_file.flush()    # write the buffer

    cmd = DNAClustCommandline(input_file=tmp_file.name,
                              similarity=identity_cutoff,
                              **kwargs)
    
    stdout, stderr = cmd()
    return DNAClustIterator(StringIO.StringIO(stdout))

def _CDHITGeneralClustHandler(sequences, identity_cutoff, cmdline, **kwargs):
    """Cluster sequences with the CD-HIT family of tools.
    
    Clusters the given sequences with the given CD-HIT command and return yield
    SeqCluster objects for each cluster. The cluster radius is given by identity_cutoff.
    Other keyword parameters are passed to the CD-HIT command.

    General function for all CD-HIT tools. Specific versions of CD-HIT should be invoked
    with CDHITClustIterator, CDHITESTClustHandler, CDHIT454ClustHandler.
    """
    # make temp. directory
    temp_dir_path = tempfile.mkdtemp()

    try:
        # write sequences to file in temp. directory
        sequence_filename = os.path.join(temp_dir_path, "seq.fasta")
        sequence_file = open(sequence_filename, "w")
        SeqIO.write(sequences, sequence_file, "fasta")
        sequence_file.close()

        cmd = cmdline(input_file=sequence_filename,
                      output_file=os.path.join(temp_dir_path, "clusters"),
                      **kwargs)
        cmd()
        
        # read all the clusters since we're going to delete everything
        clusters = list(CDHITClustIterator(open(os.path.join(temp_dir_path, "clusters.clstr"), "r")))
    finally:
        shutil.rmtree(temp_dir_path)    # remove the temp. directory

    return clusters

def CDHITESTClustHandler(sequences, identity_cutoff, **kwargs):
    """Cluster nucleotide sequences with CD-HIT-EST.
    
    Clusters the given nucleotide sequences with CD-HIT-EST and return yield SeqCluster
    objects for each cluster. The cluster radius is given by identity_cutoff. Other
    keyword parameters are passed to the CD-HIT-EST command.
    """
    # TODO check that sequences are nucleotide sequences
    if 'description_len' not in kwargs and 'd' not in kwargs:
        kwargs['description_len'] = 0
    return _CDHITGeneralClustHandler(sequences, identity_cutoff, CDHITESTCommandline, **kwargs)

def CDHIT454ClustHandler(sequences, identity_cutoff, **kwargs):
    """Cluster nucleotide sequences with CD-HIT-454.
    
    Clusters the given nucleotide sequences with CD-HIT-454 and return yield SeqCluster
    objects for each cluster. The cluster radius is given by identity_cutoff. Other
    keyword parameters are passed to the CD-HIT-454 command.
    """
    # TODO check that sequences are nucleotide sequences
    return _CDHITGeneralClustHandler(sequences, identity_cutoff, CDHIT454Commandline, **kwargs)

def CDHITClustHandler(sequences, identity_cutoff, **kwargs):
    """Cluster amino acid sequences with CD-HIT.
    
    Clusters the given amino acid sequences with CD-HIT and return yield SeqCluster
    objects for each cluster. The cluster radius is given by identity_cutoff. Other
    keyword parameters are passed to the CD-HIT command.
    """
    # TODO check that sequences are amino acid sequences
    if 'description_len' not in kwargs and 'd' not in kwargs:
        kwargs['description_len'] = 0
    return _CDHITGeneralClustHandler(sequences, identity_cutoff, CDHITCommandline, **kwargs)

if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
