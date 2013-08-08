#! /usr/bin/python

import sys
import FileUtility
import Error
import pickle

class GenomesToTaxon():
	def __init__(self):
		self.Converter={};

	def prepareConverter(self, file):
		FileUtility.isValid(file);
		sys.stderr.write("counting lines\n");
		numOfLines=FileUtility.coutLines(file);
		sys.stderr.write(str(numOfLines)+" to read\n");
		readed=0;
		for line in open(file):
			genome=int(line.split()[0]);
			taxon=int(line.split()[1]);
			self.Converter[genome]=taxon;
			readed+=1;
			if(readed%100000==0):
				sys.stderr.write(str(readed)+" lines readed on "+str(numOfLines)+"\n");
	
	def convertToTaxon(self, genome):
		if(self.genomeIsValid(genome)):
			return(self.Converter[genome]);

	def getConverter(self):
		return(self.Converter);

	def genomeIsValid(self, genome):
		detect=0;
		if(self.getConverter()[genome]):
			detect=1;
		return detect;
#to do: finish this module by dumping and loading from a pickle file.
#test performances.
	#read for a binary file with pickle
	#def fillConverter(self, file):
	#	pass;
#test main      
if __name__=="__main__":
	file=sys.argv[1];
	converter=GenomesToTaxon();
	converter.fillConverter(file);
			
		


