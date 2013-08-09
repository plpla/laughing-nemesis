#! /usr/bin/python

"""
Takes sys.argv and determine what are the parameters.
To add an option: create new variable and create the detection and getters functions.
"""

import argparse
import sys
import FileUtility

class OptionParser():
	def __init__(self, args):
		self.Parser=argparse.ArgumentParser(description="Determine the last common ancester of contigs from RayMeta.")
		
		self.subparsers=self.Parser.add_subparsers(help='sub-command help');
		#prepare
		self.parser_prepare=self.subparsers.add_parser('prepare', help='Prepare data for the construction fo the tree');
		self.parser_prepare.add_argument('-t', type=str, help='The TreeOfLife-Edges.tsv file', required=True);
		self.parser_prepare.add_argument('-n', type=str, help='The Taxon-Names.tsv file', required=True);
		self.parser_prepare.add_argument('-f', type=str, help='The GenomeToTaxon.tsv file', required=True);
		#run
		self.parser_run=self.subparsers.add_parser('run', help='Find the last common ancester');
		self.parser_run.add_argument('-d', type=str, help='The biological abundance directory', required=True);
		
		#parse the args...
		self.Arguments=vars(self.Parser.parse_args(args));

	#Needed to get the args
	def getArguments(self):
		return(self.Arguments);	
		
	def getParser(self):
		return(self.Parser);
if __name__=="__main__":
	parser=OptionParser(sys.argv[1:]);
	arg=parser.getArguments();
	print(arg);
	print(arg['f']);
	print(FileUtility.isValid(arg['f']));

#old code. Will be removed in the future.
"""
class OptionParser:
	def __init__(self, arguments):
		self.TaxonNames=detectTaxon_names(arguments);
		self.Edges=detectEdges(arguments);
		self.NumberOfThreads=detectThreads(arguments);

	def detectTaxon_names(arguments):
		taxonNames="";
		if("-t" in arguments):
			position=arguments.index("-t");
			position+=1;
			taxonNaarguments.index("-p");
			position+=1;
			threads=arguments[position];
		return taxonNames;

	def getTaxonNames(self):
		return self.TaxonNames;

	def detectEdges(arguments):
		edges=""
		if("-e" in arguments):
			position=arguments.index("-e");
			position+=1;
			edges=arguments[position];
		return edges;

	def getEdges(self):
		return self.Edges;

	def detectThreads(arguments):
		threads=1;
		if("-p" in arguments):
			position=arguments.index("-p");
			position+=1;
			threads=arguments[position];
		return threads;

	def getNumberOfThreads(self):
		return self.NumberOfThreads;
"""
