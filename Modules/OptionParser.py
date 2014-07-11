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
        self.parser = argparse.ArgumentParser(description="Set of tools to analyse metagenomic assembly produced by Ray<eta.")

        self.subparsers = self.parser.add_subparsers(help='sub-command help')

        ##########         prepare             #########################
        self.parser_prepare = self.subparsers.add_parser('prepare',
                                                         help='Prepare data for taxonomic analysis (lca)')
        self.parser_prepare.add_argument('-t', type=str, help='The TreeOfLife-Edges.tsv file', default=None,
                                         required=True)
        self.parser_prepare.add_argument('-n', type=str, help='The Taxon-Names.tsv file', default=None,
                                         required=True)
        self.parser_prepare.add_argument('-f', type=str, help='The GenomeToTaxon.tsv file', default=None,
                                         required=True)
        ##########           lca                #######################
        self.parser_lca = self.subparsers.add_parser('lca', help='Find the last common ancester of each contigs')
        self.parser_lca.add_argument('-d', type=str, help='The biological abundance directory',
                                     required=True)
        self.parser_lca.add_argument('-c', type=str, help='The contig file (.fasta)',
                                     required=True)
        self.parser_lca.add_argument('-i', type=str, help='File containing a list of contigIdentification file',
                                     required=False)
        self.parser_lca.add_argument('-b', type=int, help='Maximum number of best match to consider. ' +
                                                          'Impact compute time (default 10)',
                                     default=10, required=False)
        self.parser_lca.add_argument("-path", help="File containing the path to different ContigIdentification.tsv" +
                                     " files", type=str, required=False)
        self.parser_lca.add_argument("-r", help="Maximum taxonomic level to identify contig" +
                                                "(min = root, max= subspecies)" +
                                     "Can cause the lost of some result due to unclassified taxon",
                                     type=str, required=False)
        #The -r option will probably cause some strange behavior for 2 reasons: Unclassified taxon and taxon identified
        #at a deeper taxonomic level than asked depth

        self.parser_lca.add_argument("-v", help="Increase verbosity (dafault False)",
                                 type=bool, default=False, required=False)

        ##########        Identify                ######################

        self.parser_identification = self.subparsers.add_parser("identify",
                                                                help="Identify the best matches for each contigs")
        self.parser_identification.add_argument("-d", type=str, help='The biological abundance directory',
                                                required=True)
        self.parser_identification.add_argument("-b", type=int, help="Number of best match to keep (default 10)",
                                                default=10, required=False)
        self.parser_identification.add_argument("-path",
                                                help="File containing the path to different ContigIdentification.tsv"+
                                                     " files", type=str, required=False)


        ##########          Plot                 ##########################

        self.parser_plot = self.subparsers.add_parser("plot", help="Plot Biological abundances in a sec")

        self.parser_plot.add_argument("-f", type=str, help="A biological abundance file produced by RayMeta")
        self.parser_plot.add_argument("-m", type=float, help="Minimum taxon proportion to consider (default=0.00001",
                                      required=False, default=0.00001)
        self.parser_plot.add_argument("-t", type=str, help="File type: db or taxonomy. Dafault = taxonomy",
                                      required=False, default="taxonomy")
        self.parser_plot.add_argument("-value", type=bool, help="Show values for each bar", required=False,
                                      default=False)


        #parse the args...
        self.Arguments = vars(self.parser.parse_args(args))

    #Needed to get the args
    def getArguments(self):
        return self.Arguments

    def getParser(self):
        return self.parser

if __name__=="__main__":
    parser = OptionParser(sys.argv[1:])
    arg = parser.getArguments()
    print(arg)
    print(arg['f'])
    print(FileUtility.isValid(arg['f']))
