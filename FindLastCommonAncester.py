#! /usr/bin/python

"""
This program is intended to identify the last common ancester of the contigs.
It uses the output of MetaRay and some files that are used by Ray...
For help on how to use this program: FindLastCommonAncester.py -h
"""


import os
import Modules.OptionParser as OptionParser
import Modules.Taxon as Taxon
import Modules.GenomesToTaxon_sql as gtt
from Modules.FileUtility import *

possible_taxonomic_level=["species", "phylum", "order", "species group", "kingdom", "class", "genus", "family",
                          "subspecies"]
possible_taxonomic_level=["phylum", "class", "order", "family", "genus", "species"]


def prepare_data(args):
    """
    Construct the DB and check the tree
    """
    sys.stderr.write("Preparing data...\n")
    sys.stderr.write("Checking the tree")
    tree = Taxon.Taxonomic_tree()
    tree.read_tree_of_life(args['t'])
    tree.check_tree()
    converter = gtt.Genomes_to_taxon()
    converter.prepare_converter(args['f'], args['n'])

def prepare_tree_of_life(args):
    """
    Load the pickled tree of life
    """
    tree = Taxon.Taxonomic_tree()
    tree.read_tree_of_life(args['t'])
    return tree

def prepare_genome_taxon_converter():
    """
    Load the pickled converter
    """
    converter = gtt.Genomes_to_taxon()
    return converter

def find_contigs_ID(args):
    """"
    Similare to laughing-nemesis identify and return a structure of contigs and their identifications
    """
    sys.stderr.write("Searching the best matches for each contigs based on Ray output\n")
    if args['i']:
        contigs_id_file = readPathsFile(args['i'])
    else:
        contigs_id_file = getPathsFromDirectory(args['d'])
    checkFiles(contigs_id_file)
    contigs = readContigsTSVfile(args['d'])
    for files in contigs_id_file:
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
    if args['c']:
        sequences = read_fasta_file(args['c'])
        subset = {}
        sys.stderr.write("Before: %s" %len(contigs))
        for seq in sequences:
            subset[seq] = contigs[seq]
        contigs = subset
    sys.stderr.write("%s contigs to analyse" %len(contigs))
    return contigs

"""
Need to redesign so that special cases are sent to a file
"""
def execute_LCA(contigs, tree, converter, args):
    """
    Search the lca of contigs
    :param contigs: A dictionnary of contigs that were previously identified
    :param tree: The taxonomic tree
    :param converter: The genome to taxon converter
    :param verbosity: If true will trace the LCA search to stderr
    :return: The contigs dictionnary with the LCA.
    """
    lca = -1
    num_of_contigs = len(contigs)
    win_size = num_of_contigs / 100
    i = 0
    current = 0
    sys.stderr.write("\r%s %% contigs done" % current)
    for contig in contigs:
        #sys.stderr.write("Working on "+contig+"\n")
        contigs[contig].find_lca(tree, converter, args['v'])
        contigs[contig].clean_contig_identifications()
        contigs[contig].compute_history(converter, tree)
        if args['e'] == "valid":
            contigs[contig].get_taxonomic_lca(args['s'])
        elif args['e'] == 'total':
            contigs[contig].refine_lca(tree, converter, args['s'])
        elif args['e'] == 'level':
            if args['l'] in possible_taxonomic_level:
                contigs[contig].set_lca_by_level(args['l'])
            else:
                ValueError("The level was not specified correctly")
        else:
            raise ValueError("-e option has an invalid value")
        contigs[contig].update_LCA_name_and_rank(converter)
        i += 1
        if i % win_size == 0:
            current += 1
            sys.stderr.write("\r %s %% contigs done" % current)
    sys.stderr.write("\n")
    return contigs


#TODO:#Optimize output functions to Fred needs

def out_by_max_depth(contigs, tree, level):
    """
    Try to bring all contigs to the same taxonomic level. There are 2 possible problems:
    Taxon that have "no rank" and taxon that are deeper than the wanted level.
    For the second one, we just output their true lca level.
    :param contigs: A dict of BiologicalAbundanceContigObject
    :param tree: A full and valid TaxonomicTree. Validity is not verified
    :param level: The taxonomic level at which the user want to bring the contigs
    :return: Nothing but some string to stdout.
    """
    for contig in contigs:
        if contigs[contig].LCA_name == 'Unknown':
            print("%s\t%s\t%s\t%s\t%s" % (contig, contigs[contig].getLengthInKmer(),
                                          contigs[contig].get_coverage_detph(), contigs[contig].LCA_id,
                                          contigs[contig].LCA_name))
            continue
        taxon = tree.get_node(contigs[contig].LCA_id)
        #3 possible cases: the taxonomic level is the one selecter
        # the taxonomic level is too high (have to search for a taxon at the wanted level)
        # the taxonomic level is too low (we will find root before a taxon at the wanted level)
        t_name = converter.get_taxon_name(taxon.id)
        t_rank = converter.get_taxon_rank(taxon.id)
        if t_rank == level:
            print("%s\t%s\t%s\t%s\t%s" % (contig,contigs[contig].getLengthInKmer(),
                                          contigs[contig].get_coverage_detph(), contigs[contig].LCA_id,
                                          contigs[contig].LCA_name))
        else:
            #We search for either root or a taxon at the wanted level
            while True:
                if taxon.parent == "" or taxon.id == 1:
                    #At this point we are at the root
                    break
                taxon = tree.get_node(taxon.parent)
                t_name = converter.get_taxon_name(taxon.id)
                t_rank = converter.get_taxon_rank(taxon.id)
                if t_rank == level:
                    #At this point we have reached the wanted level
                    break
            print("%s\t%s\t%s\t%s\t%s" % (contig, contigs[contig].getLengthInKmer(),
                                          contigs[contig].get_coverage_detph(), taxon.id, t_name))

def out_with_history(contigs):
    header = "Contig-name\tContig_length_in_kmers\tContig_mode_kmer_depth\tTotal_colored_kmer\tLCA_taxon_id\tLCA_name\tLCA_rank\tLCA_score"
    for level in possible_taxonomic_level:
        header += "\t"+level+"\t"+level+"_score"
    print(header)
    for contig in contigs:
        line = ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, contigs[contig].getLengthInKmer(),
                                      contigs[contig].get_coverage_detph(), contigs[contig].coloredKmers, contigs[contig].LCA_id,
                                      contigs[contig].LCA_name, contigs[contig].LCA_rank, contigs[contig].LCA_score))
        for level in possible_taxonomic_level:
            line += "\t"+str(contigs[contig].history[level][2])+"\t"+str(contigs[contig].history[level][1])
        print(line)




def out_by_contig(contigs):
    """
    Simple tsv output 3 columns
    :param contigs: The contigs to out.
    :return: None
    """
    print("Contig-name\tContig_length_in_kmers\tContig_mode_kmer_depth\tTotal_colored_kmer\tLCA_taxon_id\tLCA_name\tLCA_rank\tLCA_score")
    for contig in contigs:
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig, contigs[contig].getLengthInKmer(),
                                      contigs[contig].get_coverage_detph(), contigs[contig].coloredKmers, contigs[contig].LCA_id,
                                      contigs[contig].LCA_name, contigs[contig].LCA_rank, contigs[contig].LCA_score))


if __name__=="__main__":
    if len(sys.argv) == 1:
        print(__doc__)
    parser = OptionParser.OptionParser(sys.argv[1:])
    args = parser.getArguments()

    if sys.argv[1] == "prepare":
        if args['t'] and args['f'] and args['n']:
            prepare_data(args)

    if sys.argv[1] == "lca" and args['d'] and args['t']:
        sys.stderr.write("Loading tree of life\n")
        tree = prepare_tree_of_life(args)
        sys.stderr.write("Tree of life loaded!\n")
        sys.stderr.write("Loading genome to taxon converter\n")
        converter = prepare_genome_taxon_converter()
        sys.stderr.write("Genome to taxon converter loaded\n")
        sys.stderr.write("Searching best matches for contigs\n")
        contigs = find_contigs_ID(args)
        sys.stderr.write("Searching is done!\n")
        sys.stderr.write("Searching LCA\n")
        contigs = execute_LCA(contigs, tree, converter, args)
        if args['o'] == "historical":
            out_with_history(contigs)
        if args['o'] == 'lca':
            out_by_contig(contigs)


