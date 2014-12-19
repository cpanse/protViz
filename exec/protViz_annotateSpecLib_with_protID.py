#!/usr/bin/python
# -*- coding: latin1 -*-
#
# Christian Panse <cp@fgcz.ethz.ch>
# Copyright (C) 2014 Functional Genomics Center Zurich ETHZ|UZH. All rights reserved.
# 
# # Licensed under  GPL version 3
#
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/exec/protViz_annotateSpecLib_with_protID.py $
# $Id: protViz_annotateSpecLib_with_protID.py 6530 2014-06-17 08:03:50Z cpanse $

import os
import sys
import re
from optparse import OptionParser

class Fasta():
    """
        reading a FASTA formated file
        the protein id and the sequence are kept
    """
    def __init__(self, fastaFilename):
        self.fastaFilename=fastaFilename

        self.proteinSeqDict=dict()
        self.proteinDescDict=dict()
        protSeqList=[]

        regexProtID=re.compile("^>.+")
        protID=''

        try:
            with open(self.fastaFilename) as f:
                for line in f:
                    # chomp
                    s=line.rstrip()
                    if (re.match(regexProtID, s)):
                        if (len(protSeqList)>0):
                            protSeqString=''.join((protSeqList))
                            if self.proteinSeqDict.has_key(protSeqString):
                                pass
                            else:
                                self.proteinSeqDict[protSeqString] = list()

                            self.proteinSeqDict[protSeqString].append(protID)

                        protSplit=s.split(" ")
                        protID=protSplit[0].replace(">", "")
                        protDesc=''.join((protSplit[1:]))
                        self.proteinDescDict[protID]=protDesc

                        protSeqList=[]
                    else:
                        protSeqList.append(s)

        except:
            print "Error: problems with reading fasta file"
            raise

class SpecLib(Fasta):
    """
    reads a file, expects a peptide sequence on the 2nd column, 
    and does a inefficient search for the sequence on the Fasta object.
    """
    def searchTrypticPeptideSeq(self, peptide):
        res=list()
        for i in self.proteinSeqDict.keys():
            if re.search("(([RK])|(^)|(^M))" + peptide, i):
                protID=self.proteinSeqDict[i]
                res.append(','.join((protID)))
        if len(res) == 0:
            return ['NA']
        return res

    def joinSpecLibWithFasta(self, speclibfilename, outputfilename):
        old_peptide=None
        try:
            with open(speclibfilename) as speclib:
                for line in speclib:
                    s=str(line).rstrip().split('\t')
                    if len(s)>1:
                        peptide=s[1]
                        if peptide == old_peptide:
                            pass
                        elif peptide == 'peptide_sequence':
                            res = ['protein_list']
                        else:
                            res=self.searchTrypticPeptideSeq(peptide)

                        with open(outputfilename, "a+") as output:
                            output.write('\t'.join((s)) + "\t" + ','.join((res)) + '\n')
                        old_peptide=peptide
        except:
            print "Error: problems with reading speclib file"
            raise

if __name__ == "__main__":
    parser = OptionParser(usage="--fasta [FASTA] --speclib [CRAN protViz::genSwathIonLib generated file]")

    parser.add_option("-f", "--fasta", dest="fastafile",
                      help="FASTA formated file", metavar="string")

    parser.add_option("-s", "--speclib", dest="speclibfile",
                      help="CRAN protViz::genSwathIonLib generated file", metavar="string")

    parser.add_option("-o", "--output", dest="outputfile",
                      help="where to put the output", metavar="string")

    (options, args) = parser.parse_args()

    if options.fastafile and  options.speclibfile and options.outputfile:
        pass
    else:
        parser.error("incorrect arguments")
        parser.print_help()
        sys.exit(1)

    sl=SpecLib(fastaFilename=options.fastafile) 
    sl.joinSpecLibWithFasta(speclibfilename=options.speclibfile, outputfilename=options.outputfile)
