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

def findLCA(args):
	sys.stderr.write("Searching the best matches for each contigs based on Ray output\n");
	
	
	
 __name__=="__main__":
	if(len(sys.argv)==1):
		print(__doc__)
	parser=OptionParser.OptionParser(sys.argv[1:]);
	args=parser.getArguments();
	if(args['t'] and args['f'] and args['n']):
		prepareData(args);
	if(args['d']):
		findLCA(args);
	

