#! /usr/bin/python

"""
Class used to construct the taxonomic tree. Use number ID.
Taxon represent a node in the tree.
"""

import threading
import logging	#for debuging
import Error
import FileUtility


class Taxon():
	def  __init__(self, taxonid):
		self.ID==taxonid;
		self.Parent="";
		self.Kids=[];
	
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

	def setParent(self, parent):
		if parent!="":
			self.Parent=parent;
"""
The taxonomic tree
"""
class TaxonomicTree():
	def __init__(self):
		self.nodes={};
	def nodeExist(self, nodeId):
		detect=0;
		try:
			self.nodes.get(nodeId);
			detect=1;
		except:
			detect=0;
		return(detect);
	def addNode(self, node):
		if(self.nodeExist(node.getId())):
			Error.error("Logic error: node already exist. Node number:"+str(node.getId()));
		self[node.getId()]=node;
	def getNode(self, nodeId):
		if(self.nodeExist(nodeId)):
			return(self.nodes[nodeId]);
		else:
			Error.error("Logic error: required node does not exist");

	def checkTree(self):
		for i in self.nodes:
			for j in self.nodes[i].getParents():
				if(self.nodeExist(j.getId())):
					continue;
				else:
					Error.error("Logic error: tree is not valid");
			for k in self.nodes[i].getKids():
				if(self.nodeExist(k.getId())):
                                        continue;
                                else:
                                        Error.error("Logic error: tree is not valid");


@staticmethod
def readTreeOfLife(file):
	FileUtility.isValid(file);
	tree=TaxonomicTree();
	for line in open(file):
		parent=int(line.split()[0]);
                kid=int(line.split()[1]);
                if(tree.nodeExist(parent)):
                        if(tree.nodeExist(kid)):
                                tree.getNode(parent).addKid(kid);
				tree.getNode(kid).setParent(parent);
                        else:   
                                newNode=Taxon(kid);
				tree.addNode(newNode);
                                tree.getNode(parent).addKid(kid);
                                tree.getNode(kid).setParent(parent);
                else:   
                        if(tree.exist(kid)):
                                newNode=Taxon(parent);
                                tree.addNode(parent);
                                tree.getNode(kid).setParent(parent);
                                tree.getNode(parent).addKid(kid);
                        else:
				newKid=Taxon(kid);
				newParent=Taxon(parent);
				newKid.setParent(parent);
				newParent.addKid(kid);
				tree.addNode(newKid);
				tree.addNode(newParent);

				
			











































	
		
		





