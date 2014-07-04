#!/usr/bin/python

from Error import *

class ContigObject(object):
    "Class used for DeNovo assembly analysis"

    def __init__(self, Name="", LengthInKmer=0, ColoredKmers=0):
        self.name=Name;
        self.lengthInKmer=LengthInKmer;
        self.coloredKmers=ColoredKmers;

    def getName(self):
        return self.name;

    def getLengthInKmer(self):
        return int(self.lengthInKmer);

    def getColoredKmers(self):
        return int(self.coloredKmers);

    def setName(self, string):
        self.name=str(string);

    def setLengthInKmers(self, length):
        self.lengthInKmer=length;

    def setColoredKmers(self, number):
        self.coloredKmers=number;

    def showTSV(self):
        line=self.name+"\t"+str(self.lengthInKmer)+"\t"+str(self.coloredKmers);
        return line;


		
