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
try:
	import cPickle as pickle
except:
	import pickle
class TaxonName():
	def __init__(self):
		self.Name="";
		self.Rank="";
	
	def getName(self):
		return(self.Name);

	def getRank(self):
		return(self.Rank);

	def setAll(self, name, rank):
		self.Name=name;
		self.Rank=rank;

	def setName(self, name):
		self.Name=name;

	def setRank(self, rank):
		self.Rank=rank;

class Taxon():
	def  __init__(self, taxonid):
		self.ID=taxonid;
		self.Parent="";
		self.Kids=[];
		self.TaxonName=TaxonName();
	
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

	def getTaxonName(self):
		return self.TaxonName;

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
	
	def dumpTree(self, file):
		sys.stderr.write("Dumping tree to disk\n");
		f=open(file, 'wb');
		pickle.dump(self.getNodes(), f, protocol=2);

	def loadTree(self, file):
		sys.stderr.write("Loading tree\n");
		f=open(file, 'rb');
		self.Nodes=pickle.load(f)
	
	#return the last common ancerstor in the tree. If it is not found, return -1 and print a warning
	def findLCA(self, node1, node2):
		if(self.nodeExist(node1) and self.nodeExist(node2)):
			pass
		else:
			Error.error("Logic error, asked to find LCA for a node that does not exists");
		parentsOfNode1=[];
		parentsOfNode1.append(node1);
		parentsOfNode2=[];
		parentsOfNode2.append(node2);
		lca=-1;
		while(self.getNodes[node1].getParents()):
			parentToAdd=self.getNodes[node1].getParents();
			parentsOfNode1.append(parentToAdd);
			node1=parentToAdd;
		while(self.getNodes[node2].getParents()):
			parentToAdd=self.getNodes[node2].getParents();
			parentsOfNode2.append(parentToAdd);
			node2=parentToAdd;
		for i in parentsOfNode1:
			if i in parentsOfNode2:
				lca=i;
				break;
		if(lca==-1):
			Errror.warning("Have not found the LCA for nodes:" +str(node1)+" and "+str(node2));
		return i
				
	def addTaxonName(self, file):
		sys.stderr.write("Adding taxon name to tree\n");
		FileUtility.isValid(file);
		for line in open(file):
			id=int(line.split('\t')[0]);
			name=line.split('\t')[1];
			rank=line.split('\t')[2];
			if(self.nodeExist(id)):
				self.getNodes()[id].getTaxonName().setAll(name, rank);
			else:
				Error.warning("A taxon is missing");
			
						
		
	def checkTree(self):
		numOfNode=len(self.getNodes());
		sys.stderr.write(str(numOfNode)+ " to be checked\n");
		checked=0;
		root=[];
		for i in self.Nodes:
			#print(self.Nodes[i].getId()==1);
			#if not(self.Nodes[i].getId()):
			#	print("my parents")
			#	print(self.Nodes[i].getParents());
			#	print("my kids")
			#	print(self.Nodes[i].getKids());
			if(self.Nodes[i].getParents()):
				if(self.nodeExist(self.Nodes[i].getParents())):
					pass;
					#print("One parent checked")					
				else:
					Error.error("Logic error: tree is not valid because this parent does not exist:"+str(self.Nodes[i].getParents()));
			else:
				sys.stderr.write("Found a root\n")
				root.append(self.Nodes[i]);
			for k in self.Nodes[i].getKids():
				if(self.nodeExist(k)):
					pass;
					#print("One kid checked");
                                else:
                                        Error.error("Logic error: tree is not valid because this kid does not exist:"+str(k));
			checked+=1;
			if(checked%50000==0):
				sys.stderr.write(str(checked)+" nodes checked\n");
		if(len(root)>1):
			Error.warning("There is more than one root\n");
		sys.stderr.write("Taxonomic tree is considered valid\n");


	def readTreeOfLife(self, file):
		FileUtility.isValid(file);
		toBeRead=FileUtility.countLines(file);
		sys.stderr.write(str(toBeRead)+" lines to read to construct tree\n");
		numOfLine=0;
		for line in open(file):
			parent=int(line.split()[0]);
			kid=int(line.split()[1]);
			if(self.nodeExist(parent)):
				if(self.nodeExist(kid)):
                                		self.getNode(parent).addKid(kid);
						self.getNode(kid).setParent(parent);
                        	else:   
                                	newNode=Taxon(kid);
					self.addNode(newNode);
                                	self.getNode(parent).addKid(kid);
                                	self.getNode(kid).setParent(parent);
                	else:   
                        	if(self.nodeExist(kid)):
                                	newNode=Taxon(parent);
                                	self.addNode(parent);
                                	self.getNode(kid).setParent(parent);
                                	self.getNode(parent).addKid(kid);
                       		else:
					newKid=Taxon(kid);
					newParent=Taxon(parent);
					newKid.setParent(parent);
					newParent.addKid(kid);
					self.addNode(newKid);
					self.addNode(newParent);
			numOfLine+=1;
			if(numOfLine%100000==0):
				sys.stderr.write(str(numOfLine)+" lines read\n");
		sys.stderr.write("Tree base constructed\n");

#test main	
if __name__=="__main__":
	files=sys.argv[1];
	sys.stderr.write("Loading taxonomic tree\n")
	tree=TaxonomicTree();
	tree.readTreeOfLife(files);
	sys.stderr.write("Taxonomic tree is loaded\n");
	sys.stderr.write("Adding taxon name to tree\n");
	tree.addTaxonName(sys.argv[2]);
	sys.stderr.write("Taxon name added to tree\n");
	sys.stderr.write("Checking if tree is valid\n");
	tree.checkTree();
	sys.stderr("Writing tree to file");
	tree.dumpTree("Data/Tree.bin");
	tree=TaxonomicTree();
	sys.stderr.write("Checking tree file");
	tree.loadTree("Data/Tree.bin");
	sys.stderr.write("Tree is ready");
	tree.checkTree();
