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

ConverterBinaryFile="Data/Converter.bin"
TreeBinaryFile="Data/Tree.bin";

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

def executeLCA(contigs):
	id1="";
	id2=""
	for id in contigs[contig].contigIdentifications:
		print(id.getSequenceName());
		
		
		
	
	

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
			contigs=findContigsID(args);
			executeLCA(contigs);

	

