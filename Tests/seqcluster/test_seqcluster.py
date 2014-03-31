#!/usr/bin/env python

from __future__ import print_function

import StringIO

from Bio.seqcluster.applications import DNAClustCommandline
from Bio.seqcluster import DNAClustIterator

from Bio.seqcluster.applications import CDHITCommandline
from Bio.seqcluster import CDHITClustIterator

cmd = DNAClustCommandline(similarity=0.8, header=True, threads=2, inputfile="test_sequences.fasta")
stdout, stderr = cmd()
clusters = DNAClustIterator(StringIO.StringIO(stdout))
for cluster in clusters:
    print(cluster.name)
    for member in cluster:
        if member == cluster.representative:
            print("\t" + member.name + "*")
        else:
            print("\t" + member.name)

print()     # blank line

cmd = CDHITCommandline(cutoff=0.8, threads=2, inputfile="test_sequences.fasta", outputfile="tmp")
stdout, stderr = cmd()
clusters = CDHITClustIterator(open("tmp.clstr", "r"))
for cluster in clusters:
    print(cluster.name)
    for member in cluster:
        if member == cluster.representative:
            print("\t" + member.name + "*")
        else:
            print("\t" + member.name)
