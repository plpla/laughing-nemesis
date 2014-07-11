__author__='pier-luc'
__date__ = '2014-07-09'
__version__ = 0.1

import sys
from Modules.OptionParser import OptionParser
from FindLastCommonAncester import *
from FindContigsIdWithBiologicalAbundance import *
from PlotBiologicalAbundance import *


if __name__ == "__main__":
    parser = OptionParser.OptionParser(sys.argv[1:])
    args = parser.getArguments()


    ########            PREPARE     #################
    if sys.argv[1] == "prepare":
        if args['t'] and args['f'] and args['n']:
            prepareData(args)


    ########           LCA              ###################
    if sys.argv[1] == "lca" and args['d'] and args['c']:
        if args['r'] is not None and not args['r'] in possible_taxonomic_level:
            sys.stderr.write("Bad taxonomic level. Possible choices are:\n %s\n" % possible_taxonomic_level)
            sys.exit(0)
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
        contigs = executeLCA(contigs, tree, converter, args["v"])
        if "r" in args and args['r'] is not None:
            out_by_max_depth(contigs, tree, args["r"])
        else:
            out_by_contig(contigs)


    ################       IDENTIFY     ###################
    if sys.argv[1] == "identify":
        directory = args["d"]
        numberOfBestMatch = args["b"]
        contigIdentificationsFiles = []
        if args["path"] is not None:
            contigIdentificationsFiles = readPathsFile(args["path"])
        else:
            contigIdentificationsFiles = getPathsFromDirectory(directory)
    #Second step: read the Contigs.tsv file from _DeNovoAssembly directory.
        biologicalAbundanceContigs = {}
        biologicalAbundanceContigs = readContigsTSVfile(directory)
    #third step: For each ContigsIdentification.tsv file, read each line and put it at the right place.
        for files in contigIdentificationsFiles:
            try:
                readContigIdentificationFiles(files, biologicalAbundanceContigs)
            except:
                warning("Unable to read: " + files)
    #Fourth step: Calculate PL-values.
        for contigs in biologicalAbundanceContigs:
            biologicalAbundanceContigs[contigs].calculatePLvalues()

    #Fifth step: Select only the top ?? identification for each contigs
        for contigs in biologicalAbundanceContigs:
            if len(biologicalAbundanceContigs[contigs].contigIdentifications) <= 0:
                continue
            else:
                biologicalAbundanceContigs[contigs].selectBestIdentifications(numberOfBestMatch)
    #Last step: Write to stdout!
        showTSV(biologicalAbundanceContigs)


    ############           PLOT_single           #########################
    if sys.argv[1] == "plot":
        data = {}
        if args["t"] == "taxonomy":
            data = read_taxonomy_file(args["f"], args["m"])
        elif args["t"] == "db":
            data = read_db_file(args["f"])
        else:
            raise ValueError("The file type specified is invalid. Must be 'taxonomy' or 'db'")
        stacked_bar_plot_single_simple(data, args["value"])



