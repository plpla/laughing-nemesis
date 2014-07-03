#! /usr/bin/python

"""
This program is intended to identify the last common ancester of the contigs.
It uses the output of MetaRay and some files that are used by Ray...
For help on how to use this program: FindLastCommonAncester.py -h
"""

import sys
import os
import Modules.OptionParser as OptionParser
import Modules.Taxon as Taxon
import Modules.GenomesToTaxon as GenomesToTaxon
from Modules.Functions import *
from Modules.FileUtility import *

######################
#Probably the only thing you could modify in the code!
ConverterBinaryFile="Data/Converter.bin"
TreeBinaryFile="Data/Tree.bin";
######################

def prepareData(args):
	sys.stderr.write("Preparing data...\n")
	tree=Taxon.TaxonomicTree();
	tree.readTreeOfLife(args['t']);
	tree.addTaxonName(args['n']);
	tree.dumpTree(TreeBinaryFile);
	tree.loadTree(TreeBinaryFile);
	tree.checkTree();
	converter=GenomesToTaxon.GenomesToTaxon();
	converter.prepareConverter(args['f']);
	converter.dumpConverter(ConverterBinaryFile);

def prepareTreeOfLife():
	tree=Taxon.TaxonomicTree();
	tree.loadTree(TreeBinaryFile);
	return tree;

def prepareGenomeToTaxonConverter():
	converter=GenomesToTaxon.GenomesToTaxon();
        converter.loadConverter(ConverterBinaryFile);
	return converter
	
def findContigsID(args):
	sys.stderr.write("Searching the best matches for each contigs based on Ray output\n");
	if(args['i']):
		contigsIDfile=readPathsFile(args['i']);
	else:
		contigsIDfile=getPathsFromDirectory(args['d']);
	checkFiles(contigsIDfile);
	contigs=readContigsTSVfile(args['d']);
	for files in contigsIDfile:
		try:
			readContigIdentificationFiles(files, contigs);
		except:
			warning("Unable to read: "+ files);
	for contig in contigs:
		contigs[contig].calculatePLvalues();
	for contig in contigs:
		if(len(contigs[contig].contigIdentifications)<=0):
			continue;
		else:
			contigs[contig].selectBestIdentifications(args['b']);
	return(contigs);

"""
Need to redesign so that special cases are sent to a file
"""
def executeLCA(contigs, tree, converter):
	lca=-1;
	for contig in contigs:
		idList=contigs[contig].contigIdentifications;
		numberOfId=len(idList);
		#3 cases: 0 id, 1 id and 2 id or more
		if(numberOfId==0):
			pass;	#TODO: Will have to do something with uncoloried contigs. Are they still in there?
		if(numberOfId==1):
			sys.stderr.write("Case where there is only 1 match:\n");
			if(idList[0].getSequenceName().split('|')[0]=="gi"):
				sys.stderr.write(idList[0].getSequenceName());
				id=int(idList[0].getSequenceName().split('|')[1]);
				if(converter.genomeIsValid(id)):
					contigs[contig].LCA_id=converter.convertToTaxon(id);
				##	sys.stderr.write("Converted");
					node=tree.getNode(contigs[contig].LCA_id);
				#	sys.stderr.write("Found node");
					contigs[contig].LCA_name=node.getTaxonName().getName();
				else:
					sys.stderr.write("NOT VALID\n")
			else:
				sys.stderr.write(contigs[contig].getName()); #DERNIERE LIGNE MODIFIE
		if(numberOfId>1):
			sys.stderr.write("Case where there is more than one match\n");
			index=0;
			lca=0;
			while(lca==0):
				try:
					sys.stderr.write("Searching a lca start with genome:\n");
					sys.stderr.write(idList[index].getSequenceName()+"\n");
					lca=int(idList[index].getSequenceName().split('|')[1]);
					if not(converter.genomeIsValid(lca)):
						lca=0;
				except:
					index+=1;
					lca=0;
			#endWhile
			#LCA has now a valid value. We can iterate on each entry to find the true lca!
			#Got a problem in this section with taxon/genome id.
			id1=converter.convertToTaxon(lca);
			for identification in idList:
				if(identification.getSequenceName().split('|')[0]=="gi"):
					sys.stderr.write("In a true LCA case. Converting an id");
					id2=int(identification.getSequenceName().split('|')[1]);
					sys.stderr.write("In a true LCA case. Checking if valid");
					if(converter.genomeIsValid(id2)):
						sys.stderr.write("Searching the LCA in the tree\n");
						id2=converter.convertToTaxon(id2);
						lca=tree.findLCA(id1, id2);
						id1=lca;
				else:
					pass;#We should do something about it...
			print(lca);
				
		#id1="";
		#id2="";


if __name__=="__main__":
	if(len(sys.argv)==1):
		print(__doc__)
	parser=OptionParser.OptionParser(sys.argv[1:]);
	args=parser.getArguments();
	if(sys.argv[1]=="prepare"):
		if(args['t'] and args['f'] and args['n']):
			prepareData(args);
	if(sys.argv[1]=="run"):
		if(args['d'] and args['c']):
			sys.stderr.write("Loading tree of life\n");
			tree=prepareTreeOfLife();
			sys.stderr.write("Tree of life loaded!\n");
			sys.stderr.write("Loading genome to taxon converter\n");
			converter=prepareGenomeToTaxonConverter();
		
			sys.stderr.write("Genome to taaxon converter loaded\n");
			sys.stderr.write("Searching best matches for contigs\n");
			contigs=findContigsID(args);
			sys.stderr.write("Searching is done!\n");
			sys.stderr.write("Searching LCA\n");
			executeLCA(contigs, tree, converter);

	

