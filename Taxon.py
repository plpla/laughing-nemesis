#! /usr/bin/python

"""
Class used to construct the taxonomic tree. Use number ID.
"""

class Taxon:
	def __init__(self, taxonid, parent="", kid=""):
		self.ID==taxonid;
		self.Parent="";
		self.Kids=[];
		if(kid!=""): #not sure if its necessary to test that in here.
			self.Kids.append(kid);
	
	def hasParent(self, taxonId):
		a=0;
		if(taxonId in self.Parent):
			a=1;
		return a;

	def hasKid(sef, taxonId):
		a=0;
		if(taxonId in self.Kids):
			a=1;
		return a;

	def getKids(self):
		return self.Kids;

	def getParents(self):
		return self.parents;

	def getId(self):
		return self.ID;

	def addKid(self, taxonId):
		if not (self.hasKid(taxonId)):
			self.Kids.append(taxonId);
		else:
			print(taxonId+" is already a kid");






