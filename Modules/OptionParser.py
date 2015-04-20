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

        self.parser_lca.add_argument('-c', type=str, help='The contig file (.fasta) to search only for these contigs (not implemented yet)',
                                     required=False)

        self.parser_lca.add_argument('-i', type=str, help='File containing a list of contigIdentification file. This could produce strange results',
                                     required=False)

        self.parser_lca.add_argument('-t', type=str, help='The TreeOfLife-Edges.tsv file. Usually located in the Taxonomy repository.',
                                     default=None, required=True)

        self.parser_lca.add_argument('-b', type=int, help='Maximum number of best match to consider. ' +
                                                          'Impact compute time. (default 1000)',
                                     default=1000, required=False)

        self.parser_lca.add_argument('-o', type=str, help='Stdout format. NOT stdout file name. Default: historical',
                                     default="historical", choices=['lca', 'historical'], required=False)

        self.parser_lca.add_argument('-e', type=str, help="Type of LCA to return. valid, total or level. Level is required for -l. Default: valid",
                                    default='valid', required=False, choices=['valid', 'total', 'level'])

        self.parser_lca.add_argument("-l", help="Taxonomic level at which the LCA will be called." +
                                                "(min = phylum, max= species)." +
                                     "Can cause the lost of some result due to unclassified taxon. " +
                                     "Will return the match with the highest scrore if match at level is higher than the LCA",
                                     type=str, required=False,
                                     choices=["phylum", "class", "order", "family", "genus", "species"])

        self.parser_lca.add_argument("-s",
                                     help="Minimal score [0-1] for LCA call. See developper for info on score (will eventually be in the doc). Default: 0.9",
                                     type=float, required=False, default= 0.9)
        #The -s option will probably cause some strange behavior for 2 reasons: Unclassified taxon and taxon identified
        #at a deeper taxonomic level than asked depth

        self.parser_lca.add_argument("-v", help="Increase verbosity. Not fully implemented. Default: False",
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
        self.parser_plot.add_argument("-multi", type=bool, help="If true, the -f file contains a list of biological"+
                                                                "abundance file to plot in the same figure",
                                      required=False, default=False)
        self.parser_plot.add_argument("-stacked", type=bool, help="Create a stacked bar graph instead of the"+
                                                                  "normal graph", required=False, default=False)
        self.parser_plot.add_argument("-o", type=str, help="Name of the output file. Will not open a new window",
                                      required=False, default=None)


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
