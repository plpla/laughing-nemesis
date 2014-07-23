#!/usr/bin/python

class ContigIdentificationObject(object):
    "Contains some information about a contig that helps to its identication"
    def __init__(self, SequenceName="", SequenceLengthInKmer=0, MatchesInContig=0):
        self.sequenceName = SequenceName
        self.sequenceLenghtInKmer = SequenceLengthInKmer
        self.matchesInContig = MatchesInContig
        self.PLvalue = 0

    def setSequenceName(self, name):
        self.sequenceName = name

    def setSequenceLengthInKmer(self, length):
        self.sequenceLenghtInKmer = length

    def setMatchesInContig(self, number):
        self.matchesInContig = number

    def setPLvalue(self, value):
        self.PLvalue = value

    def getSequenceName(self):
        return self.sequenceName

    def getSequenceLengthInKmer(self):
        return self.sequenceLenghtInKmer

    def getMatchesInContig(self):
        return self.matchesInContig

    def getPLvalue(self):
        return self.PLvalue

    def showTSV(self):
        print(self.sequenceName+"\t"+self.sequenceLenghtInKmer+"\t"+self.matchesInContig)

    def calculatePLvalue(self, LengthInKmer, ColoredKmers):
        if self.matchesInContig>0 and ColoredKmers>0 and LengthInKmer>0:
            #print("using real pl-value");
            #print("MatchesInContig"+str(self.matchesInContig)+"ColoredKmer"+str(ColoredKmers)+"length"+str(LengthInKmer))
            plvalue = ((self.matchesInContig*self.matchesInContig)/ColoredKmers)/float(LengthInKmer)
        #print(plvalue);
        else:
            plvalue=0
        self.setPLvalue(plvalue)
