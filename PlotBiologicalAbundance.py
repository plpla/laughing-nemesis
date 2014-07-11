__author__='pier-luc'
__date__ = "2014-07-09"

import sys
from Modules.OptionParser import OptionParser
import numpy as np
import matplotlib.pyplot as mpl




#mpl.use('GTK')


def read_taxonomy_file(file_name, minimum_proportion):
    data = {}
    for line in open(file_name):
        if line[0] == '#':
            continue
        taxon_name = line.split("\t")[1]
        taxon_proportion = convert_exponential_notation(line.split("\t")[3])
        if taxon_proportion >= minimum_proportion:
            data[taxon_name] = taxon_proportion
    return data

def convert_exponential_notation(exponential_exp):
    if "e" in exponential_exp:
        num = float(exponential_exp.split("e")[0])
        expo = int(exponential_exp.split('e')[1])
        converted_number = num*(10**expo)
    else:
        converted_number = float(exponential_exp)
    return converted_number


def stacked_bar_plot_single_simple(data, show_value):
    number_of_sample = len(data)
    taxon_name = []
    taxon_value = []
    for taxon in data:
        taxon_name.append(taxon)
        taxon_value.append(data[taxon])
    ind = np.arange(number_of_sample)
    width = 0.75    # could be optimized
    plot = mpl.bar(ind, taxon_value, width, color='g')
    mpl.ylabel("Proportion")
    mpl.title("Taxonomic repartition")
    locs, labels = mpl.xticks(ind+width/2., taxon_name)
    mpl.setp(labels, rotation=90)
    mpl.yticks(np.arange(0, max(taxon_value), max(taxon_value)/10))
    mpl.tight_layout()
    if show_value:
        for pos, value in zip(ind, taxon_value):
            mpl.text(pos + 0.5*width, value, str(value)[0:6], ha='center', va='bottom', rotation=90)
    mpl.show()




        #for axis rotation:
#http://stackoverflow.com/questions/10998621/rotate-axis-text-in-python-matplotlib


def read_db_file(file_name):
    raise NotImplementedError


if __name__=="__main__":
    """
    New tool to plot Biological abundance directly from Ray output.
    Objective is to obtain a graph in less than a second
    Futur objective, a web page with possibly some parameters for sample exploration
    """
    parser = OptionParser(sys.argv[1:])
    args = parser.getArguments()
    args = parser.getArguments()
    if sys.argv[1] != "plot":
        raise RuntimeError("The only option available in this script is 'plot'")
    data = {}
    if args["t"] == "taxonomy":
        data = read_taxonomy_file(args["f"], args["m"])
    elif args["t"] == "db":
        data = read_db_file(args["f"])
    else:
        raise ValueError("The file type specified is invalid. Must be 'taxonomy' or 'db'")
    stacked_bar_plot_single_simple(data, args["value"])






