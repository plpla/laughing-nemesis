#!/usr/bin/python

from ContigObject import *
from ContigIdentificationObject import *

class BiologicalAbundanceContigObject(ContigObject):
	"Class representing a contig and the information about it from RayMeta"
	def __init__(self, name="", lengthInKmer=0, coloredKmers=0):
		ContigObject.__init__(self, name, lengthInKmer, coloredKmers);
		self.contigIdentifications=[];
	
	def addNewContigIdentification(self, Name, Length, Matches):
		newContigIdentification=ContigIdentificationObject(Name, Length, Matches);
		self.contigIdentifications.append(newContigIdentification);
		#print("addNewContigIdentification:"+str(len(self.contigIdentifications)));

	def calculatePLvalues(self):
		for entry in self.contigIdentifications:
			entry.calculatePLvalue(ContigObject.getLengthInKmer(self), ContigObject.getColoredKmers(self));

	def selectBestIdentifications(self, numberOfBestId):
		#New step: removePL-value of 0;
		nonNulContigIdentification=[];
		#print("Select before non-nul");
		#print(len(self.contigIdentifications));
		for entry in self.contigIdentifications:
			entry.calculatePLvalue(ContigObject.getLengthInKmer(self), ContigObject.getColoredKmers(self));
			#print(entry.getPLvalue());
			#print(ContigObject.getLengthInKmer(self));
			#print(ContigObject.getColoredKmers(self));
			if(entry.getPLvalue()>0):
				nonNulContigIdentification.append(entry);
		self.contigIdentifications=nonNulContigIdentification;
		#print("Select after select non-=null");
                #print(len(self.contigIdentifications));
		if(len(self.contigIdentifications)<=0):
			return;
		else:
			(self.contigIdentifications).sort(key=lambda x: x.PLvalue, reverse=True);
			#print("Select after sort");
                	#print(len(self.contigIdentifications));
			if(len(self.contigIdentifications)>numberOfBestId):
				self.contigIdentifications=self.contigIdentifications[0:numberOfBestId];

	def removeNulIdentification(self):
		contigIdModified=[];
		for entry in self.contigIdentifications:
			if entry.getPLvalue()>0:
				contigIdModified.append(entry);
		return contigIdModified;

#	def removeHumanIdentifications(self):
#		for entry in self.contigIdentifications:
#			if("Homo sapiens chromosome" in entry.getSequenceName()):
			
		

	def showTSV(self):
		#self.calculatePLvalues();
		if(len(self.contigIdentifications)!=0):
			#print("Longueur contigIdentifications:"+str(len(self.contigIdentifications)));
			for entry in self.contigIdentifications:
				line=ContigObject.showTSV(self);
				line=line+"\t"+entry.getSequenceName()+"\t"+str(entry.getSequenceLengthInKmer())+"\t"+str(entry.getMatchesInContig())+"\t"+str(entry.getPLvalue());
				print(line);

