#!/usr/bin/python

from Error import *

class ContigObject(object):
    "Class used for DeNovo assembly analysis"

    def __init__(self, name="", length_in_kmer=0, colored_kmers=0, coverage_depth=0):
        self.name = name
        self.lengthInKmer= length_in_kmer
        self.coloredKmers= colored_kmers
        self.mode_kmer_coverage_depth = coverage_depth

    def getName(self):
        return self.name

    def getLengthInKmer(self):
        return int(self.lengthInKmer)

    def getColoredKmers(self):
        return int(self.coloredKmers)

    def setName(self, string):
        self.name = str(string)

    def setLengthInKmers(self, length):
        self.lengthInKmer = length

    def setColoredKmers(self, number):
        self.coloredKmers = number

    def get_coverage_detph(self):
        return int(self.mode_kmer_coverage_depth)

    def showTSV(self):
        line = self.name+"\t"+str(self.lengthInKmer)+"\t"+str(self.coloredKmers)
        return line


		
