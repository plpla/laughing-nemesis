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
TreeBinaryFile = "Data/Tree.bin"
######################

def prepareData(args):
    sys.stderr.write("Preparing data...\n")
    tree = Taxon.TaxonomicTree()
    tree.readTreeOfLife(args['t'])
    tree.addTaxonName(args['n'])
    tree.dumpTree(TreeBinaryFile)
    tree.loadTree(TreeBinaryFile)
    tree.checkTree()
    converter = GenomesToTaxon.GenomesToTaxon()
    converter.prepareConverter(args['f'])
    converter.dumpConverter(ConverterBinaryFile)

def prepareTreeOfLife():
    tree = Taxon.TaxonomicTree()
    tree.loadTree(TreeBinaryFile)
    return tree

def prepareGenomeToTaxonConverter():
    converter = GenomesToTaxon.GenomesToTaxon()
    converter.loadConverter(ConverterBinaryFile)
    return converter

def findContigsID(args):
    sys.stderr.write("Searching the best matches for each contigs based on Ray output\n")
    if args['i']:
        contigsIDfile = readPathsFile(args['i'])
    else:
        contigsIDfile = getPathsFromDirectory(args['d'])
    checkFiles(contigsIDfile)
    contigs = readContigsTSVfile(args['d'])
    for files in contigsIDfile:
        try:
            readContigIdentificationFiles(files, contigs)
        except:
            warning("Unable to read: " + files)
    for contig in contigs:
        contigs[contig].calculatePLvalues()
    for contig in contigs:
        if len(contigs[contig].contigIdentifications) <= 0:
            continue
        else:
            contigs[contig].selectBestIdentifications(args['b'])
    return contigs

"""
Need to redesign so that special cases are sent to a file
"""
def executeLCA(contigs, tree, converter, verbosity):
    lca = -1
    for contig in contigs:
        if verbosity:
            sys.stderr.write("Computing %s\n" %contig)
        idList = contigs[contig].contigIdentifications
        numberOfId = len(idList)
        #3 cases: 0 id, 1 id and 2 id or more
        if numberOfId == 0:
            contigs[contig].LCA_id = "No color"
            contigs[contig].LCA_name = "Unknown"
        if numberOfId == 1:
            if verbosity:
                sys.stderr.write("Case where there is only 1 match:\n")
            if idList[0].getSequenceName().split('|')[0] == "gi":
                if verbosity:
                    sys.stderr.write(idList[0].getSequenceName()+"\n")
                id = int(idList[0].getSequenceName().split('|')[1])
                if converter.genomeIsValid(id):
                    contigs[contig].LCA_id = converter.convertToTaxon(id)
                    node = tree.getNode(contigs[contig].LCA_id)
                    contigs[contig].LCA_name = node.getTaxonName().getName()
                else:
                    contigs[contig].LCA_id = "Invalid id"
                    contigs[contig].LCA_name = "Unknown"
                    if verbosity:
                        sys.stderr.write("NOT VALID\n")
            else:
                contigs[contig].LCA_id = "No gi"
                contigs[contig].LCA_name = "Unknown"
                sys.stderr.write(contigs[contig].getName())    # DERNIERE LIGNE MODIFIE
        if numberOfId > 1:
            if verbosity:
                sys.stderr.write("Case where there is %s match\n" % numberOfId)
            index = 0
            lca = 0
            #1- We set lca to a value that is know to exist in the tree
            while lca == 0 and index < numberOfId:
                if verbosity:
                    sys.stderr.write("Searching LCA starting point with genome: ")
                    sys.stderr.write(idList[index].getSequenceName()+"\n")
                try:
                    lca = int(idList[index].getSequenceName().split('|')[1])
                except IndexError:
                    lca = 0
                    index += 1
                    continue    #TODO: This section is messy... Do some clean code bro!
                validity = converter.genomeIsValid(lca)
                if verbosity:
                    sys.stderr.write("id %s is valid: %s\n" % (lca, validity))
                if not validity:
                    lca = 0
                    index += 1
            if not converter.genomeIsValid(lca):
                contigs[contig].LCA_id = "Unidentified"
                contigs[contig].LCA_name = "Unknown"
                # If you are here, you have reached the last possible identification and
                #  have not found one that is valid.
            else:
            #LCA has now a valid value. We can iterate on each entry to find the true lca!
            #Got a problem in this section with taxon/genome id.
                id1 = converter.convertToTaxon(lca)
                for identification in idList:
                    if identification.getSequenceName().split('|')[0]=="gi":
                        if verbosity:
                            sys.stderr.write("Found a GI. Fetching name\n")
                        id2 = int(identification.getSequenceName().split('|')[1])
                        if verbosity:
                            sys.stderr.write("Converting name to taxon id. Checking if valid\n")
                        if converter.genomeIsValid(id2):
                            id2 = converter.convertToTaxon(id2)
                            if verbosity:
                                sys.stderr.write("Searching the LCA in the tree for %s and %s\n" % (id1, id2))
                            lca = tree.findLCA(id1, id2)
                            if verbosity:
                                sys.stderr.write("Resulting LCA is %s\n" %lca)
                            id1 = lca
                        else:
                            if verbosity:
                                sys.stderr.write("Sequence %s cant be converted\n" % id2)
                            continue
                    else:
                        if verbosity:
                            sys.stderr.write("%s , %s there is no GI \n" % (contig, identification.getSequenceName()))
                        continue #We should do something about it...
            contigs[contig].LCA_id = lca
            if tree.nodeExist(contigs[contig].LCA_id):
                node = tree.getNode(contigs[contig].LCA_id)
                if verbosity:
                    sys.stderr.write("Setting %s LCA_name to %s by taxon %s\n" % (contig, node.getTaxonName().getName(),
                                                                                  node.ID))
                contigs[contig].LCA_name = node.getTaxonName().getName()
            else:
                contigs[contig].LCA_name = "Unknown"
                #contigs[contig].LCA_name = name
        print("%s\t%s\t%s" % (contig, contigs[contig].LCA_id, contigs[contig].LCA_name))

        #id1="";
        #id2="";


if __name__=="__main__":
    if len(sys.argv) == 1:
        print(__doc__)
    parser = OptionParser.OptionParser(sys.argv[1:])
    args = parser.getArguments()
    if sys.argv[1] == "prepare":
        if args['t'] and args['f'] and args['n']:
            prepareData(args)
    if sys.argv[1] == "lca" and args['d'] and args['c']:
        sys.stderr.write("Loading tree of life\n")
        tree = prepareTreeOfLife()
        sys.stderr.write("Tree of life loaded!\n")
        sys.stderr.write("Loading genome to taxon converter\n")
        converter = prepareGenomeToTaxonConverter()
        sys.stderr.write("Genome to taaxon converter loaded\n")
        sys.stderr.write("Searching best matches for contigs\n")
        contigs = findContigsID(args)
        sys.stderr.write("Searching is done!\n")
        sys.stderr.write("Searching LCA\n")
        executeLCA(contigs, tree, converter, args["v"])

	

