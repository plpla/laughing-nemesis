#!/usr/bin/python

"""
This script determines the X best matches from the DataBase used with RayMeta for each Contigs of an assembly.
Usage:
python FindContigsIdWithBiologicalAbundance.py BiologicalAbundanceDirectory numberOfBestMatch FileContainingThePathToContigsIdentification.tsvFiles(optional)
If a file with the directory of ContigsIdentification.tsv files  is given, only those files will be considered.
The program display only entries for wich there is a PL-value higher than 0. It can be disabled by modifying a fct in Biological....py
"""

from Modules.BiologicalAbundanceContigObject import *
from Modules.FileUtility import *
def showTSV(biologicalAbundanceContigs):
        header = "Contig name"+'\t'+"Contig length in k-mer"+'\t'+"Contig colored k-mer"+'\t'+\
               "Sequence name"+'\t'+"Sequence length in k-mer"+"\t"+"Matches in sequence"+'\t'+"pl-value"
        print(header)
        for contigs in biologicalAbundanceContigs:
                biologicalAbundanceContigs[contigs].showTSV()



if __name__ == "__main__":
    testValidity()
    directory = sys.argv[1]
    numberOfBestMatch = int(sys.argv[2])
    if len(sys.argv) == 3:
        pathToContigsIdFile = ""
    else:
        pathToContigsIdFile = sys.argv[3]
    #First step: read the Path file.
    contigIdentificationsFiles = []
    if pathToContigsIdFile != "":
        contigIdentificationsFiles = readPathsFile(pathToContigsIdFile)
    else:
        contigIdentificationsFiles=getPathsFromDirectory(directory)
    #Second step: read the Contigs.tsv file from _DeNovoAssembly directory.
    biologicalAbundanceContigs = {}
    biologicalAbundanceContigs = readContigsTSVfile(directory)
    #third step: For each ContigsIdentification.tsv file, read each line and put it at the right place.
    for files in contigIdentificationsFiles:
        try:
            readContigIdentificationFiles(files, biologicalAbundanceContigs)
        except:
            raise IOError("Unable to read: " + files)
    #Fourth step: Calculate PL-values.
    for contigs in biologicalAbundanceContigs:
        biologicalAbundanceContigs[contigs].calculatePLvalues()

    #Fifth step: Select only the top ?? identification for each contigs
    for contigs in biologicalAbundanceContigs:
        if len(biologicalAbundanceContigs[contigs].contigIdentifications) <= 0:
            continue
        else:
            biologicalAbundanceContigs[contigs].selectBestIdentifications(numberOfBestMatch)
    #Last step: Write to stdout!
    showTSV(biologicalAbundanceContigs)
