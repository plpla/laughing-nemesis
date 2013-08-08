#! /usr/bin/python

"""
Class used to construct the taxonomic tree. Use number ID.
Taxon represent a node in the tree.
"""

import threading
import logging	#for debuging
import Error
import FileUtility
import sys

class Taxon():
	def  __init__(self, taxonid):
		self.ID=taxonid;
		self.Parent="";
		self.Kids=[];
	
	def hasParent(self, taxonId):
		a=0;
		if(taxonId in self.Parent):
			a=1;
		return a;

	def hasKid(self, taxonId):
		a=0;
		if(not self.Kids):
			if(taxonId in self.Kids):
				a=1;
		return a;

	def getKids(self):
		return self.Kids;

	def getParents(self):
		return self.Parent;

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
		self.Nodes={};
	
	def getNodes(self):
		return(self.Nodes);

	def nodeExist(self, nodeId):
		detect=0;
		try:
			self.Nodes[nodeId];
			detect=1;
		except:
			detect=0;
		return(detect);

	def addNode(self, node):
		if(self.nodeExist(node.getId())):
			Error.error("Logic error: node already exist. Node number:"+str(node.getId()));
		self.Nodes[node.getId()]=node;

	def getNode(self, nodeId):
		if(self.nodeExist(nodeId)):
			return self.Nodes[nodeId];
		else:
			Error.error("Logic error: required node does not exist");

	def checkTree(self):
		numOfNode=len(self.getNodes());
		sys.stderr.write(str(numOfNode)+ " to be checked\n");
		checked=0;
		root=[];
		for i in self.Nodes:
			#print(self.Nodes[i].getId());
			if(self.Nodes[i].getId()==1):
				print("my parents")
				print(self.Nodes[i].getParents());
				print("my kids")
				print(self.Nodes[i].getKids());
			if(self.Nodes[i].getParents()):
				if(self.nodeExist(self.Nodes[i].getParents())):
					pass;
					#print("One parent checked")					
				else:
					Error.error("Logic error: tree is not valid because this parent does not exist:"+str(self.Nodes[i].getParents()));
			else:
				sys.stderr.write("Found a node without parents\n")
				root.append(self.Nodes[i].getParents());
			for k in self.Nodes[i].getKids():
				if(self.nodeExist(k)):
					pass;
					#print("One kid checked");
                                else:
                                        Error.error("Logic error: tree is not valid because this kid does not exist:"+str(k));
			checked+=1;
			if(checked%1000==0):
				sys.stderr.write(str(checked)+" nodes checked\n");
		for b in root:
			sys.stderr.write(str(b)+"\n");
		sys.stderr.write("Taxonomic tree is considered valid\n");


def readTreeOfLife(file):
	FileUtility.isValid(file);
	tree=TaxonomicTree();
	numOfLine=0;
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
                        if(tree.nodeExist(kid)):
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
		numOfLine+=1;
		if(numOfLine%100000==0):
			sys.stderr.write(str(numOfLine)+" lines read\n");
	return(tree);

#test main	
if __name__=="__main__":
	files=sys.argv[1];
	sys.stderr.write("Loading taxonomic tree\n")
	tree=readTreeOfLife(files);
	sys.stderr.write("Taxonomic tree is loaded\n");
	sys.stderr.write("Checking if tree is valid\n")
	tree.checkTree();
			
