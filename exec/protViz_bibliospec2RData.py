#!/usr/bin/python
# -*- coding: latin1 -*-


# Christian Panse <cp@fgcz.ethz.ch>

# Copyright (C) 2014 Functional Genomics Center Zurich ETHZ|UZH. All rights reserved.
# 
# # Licensed under  GPL version 3
# $HeadURL: http://fgcz-svn.unizh.ch/repos/fgcz/testing/proteomics/R/protViz/exec/protViz_bibliospec2RData.py $
# $Id: protViz_bibliospec2RData.py 6607 2014-08-13 08:04:27Z cpanse $


"""
# INPUT: bibliospec sqlite files
# OUTPUT: extracts MS data and makes it R readable
"""

import sqlite3
import zlib
import sys
import numpy
import os
import math
import re
import string

def bibliospec2R(sqllitedb):
    count=0
    myquery = "SELECT numPeaks, peakMZ, peakIntensity, peptideSeq, precursorCharge, precursorMZ, retentionTime, peptideModSeq, score, SpectrumSourceFiles.fileName FROM SpectrumSourceFiles, RefSpectraPeaks, RefSpectra WHERE RefSpectra.id=RefSpectraPeaks.RefSpectraID and SpectrumSourceFiles.id = RefSpectra.fileID;"
    myconnect = sqlite3.connect(sqlitedb)
    mycursor = myconnect.cursor()

    myRDataName = os.path.basename(sqllitedb)

    myRFile = open(sqllitedb+".R", 'w')

    myRFile.write(myRDataName + "<-list()")

    for row in mycursor.execute(myquery):
        count = count + 1

        try:
            # COMPRESSED 
            myblob = zlib.decompress(row[1]) 
            MZ = numpy.fromstring(myblob, dtype='d')
        except zlib.error as e:
            #  NOT COMPRESSED 
            MZ = numpy.fromstring(row[1], dtype='d')
        except:
            print "+++ unexpected error +++"
            sys.exit(1)

        try:
            # COMPRESSED 
            myblob = zlib.decompress(row[2]) 
            Intensity=numpy.fromstring(myblob, dtype='f')
        except zlib.error as e:
            # NOT COMPRESSED 
            Intensity=numpy.fromstring(row[2], dtype='f')
        except:
            print "+++ unexpected error +++"
            sys.exit(1)

        myRFile.write( "\n\n\n" + myRDataName + "[[" + '{}'.format(count) + "]] <- list(" )

        myRFile.write("\n\tmZ=c(")
        for idx in range(len(MZ)):
            myRFile.write('{}'.format(MZ[idx]))
            if idx < len(MZ)-1:
                myRFile.write(", ")
        myRFile.write("),\n")

        myRFile.write("\n\tintensity=c(")
        for idx in range(len(MZ)):
            myRFile.write('{}'.format(Intensity[idx]))
            if idx < len(MZ)-1:
                myRFile.write(", ")
        myRFile.write("),\n")

        myRFile.write("\tpeptideSequence='" + str(row[3]) + "',\n")
        myRFile.write("\tcharge=" + str(row[4]) + ",\n")
        myRFile.write("\tpepmass=" + str(row[5]) + ",\n")

        myRFile.write("\tpeptideModSeq='"  + str(row[7]) + "',\n")

        varModification = row[7]
        varModification = re.sub("[A-Z]", "0,", varModification)
        varModification = re.sub("0,\[", "", varModification)
        varModification = re.sub("\]", ",", varModification)
        varModification = re.sub(",$", "", varModification)
        varModification = "c(" + varModification  + ")"
        myRFile.write("\tvarModification="  + str(varModification) + ",\n")

        # there is no mascot score defined. it looks more like an E-value
        myRFile.write("\tmascotScore="  + str( round(-10 * math.log(1E-6 + ((float(row[8])))),2)) + ",\n")

        myRFile.write("\tproteinInformation=''"  + ",\n")

        fileName = str(row[9])
        myRFile.write("\tfileName=" + repr(fileName) + ",\n")

        myRFile.write("\trtinseconds=" + str(60 * row[6]) + "\n")

        myRFile.write("); \n \n")

    # todo run R 
    myRFile.write("save(" + myRDataName + ", file='" + myRDataName + ".RData', compress=TRUE)")

    print ("run:\n R --no-save < " + sqllitedb + ".R" + "\n to geneated a '" + myRDataName + ".RData'")

def usage():
    print ("usage:")
    print ("protViz_bibliospec2RData.py file.blib")
    sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) is not 2:
        usage()

    sqlitedb=sys.argv[1]

    if os.path.isfile(sqlitedb) and os.access(sqlitedb, os.R_OK):
        bibliospec2R(sqlitedb)
    else: 
        print ("file not found or not accessable. abord.")
        sys.exit(1)

        
    print "System exit 0"

#    try:
#        t2=numpy.fromstring(row[2], dtype='float32')
#        print "size = ", sys.getsizeof(t2)
#    except:
#        pass

