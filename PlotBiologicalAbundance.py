__author__='pier-luc'
__date__ = "2014-07-09"

import sys
from Modules.OptionParser import OptionParser
import numpy as np
import matplotlib as mpl
from matplotlib.font_manager import FontProperties



#mpl.use('GTK')


def read_taxonomy_file(file_name, minimum_proportion):
    """
    Read a 4 column file produced by Ray.
    Columns are: #TaxonIdentifier        TaxonName       TaxonRank       TaxonProportion
    :param file_name: The name of the file to load (str)
    :param minimum_proportion: Minimum proportion required to plot a taxon
    :return: A dictionnary of string and values
    """
    data = {}
    for line in open(file_name):
        if line[0] == '#':
            continue
        taxon_name = line.split("\t")[1]
        taxon_proportion = convert_exponential_notation(line.split("\t")[3])
        if taxon_proportion >= minimum_proportion:
            data[taxon_name] = taxon_proportion
    return data


def read_db_file(file_name, minimum_proportion):
    """
    Read a 2 column file produced by Ray.
    Columns are: #Name   Proportion
    :param file_name: The name of the file to load (str)
    :param minimum_proportion: Minimum proportion required to plot a taxon
    :return: A dictionnary of string and values
    """
    data = {}
    for line in open(file_name):
        if line[0] == '#':
            continue
        taxon_name = line.split("\t")[0]
        taxon_proportion = convert_exponential_notation(line.split("\t")[1])
        if taxon_proportion >= minimum_proportion:
            data[taxon_name] = taxon_proportion
    return data



def convert_exponential_notation(exponential_exp):
    """
    Convert a string representing a float to a float
    :param exponential_exp: Sting of numbers and 'e'
    :return:A float
    """
    if "e" in exponential_exp:
        num = float(exponential_exp.split("e")[0])
        expo = int(exponential_exp.split('e')[1])
        converted_number = num*(10**expo)
    else:
        converted_number = float(exponential_exp)
    return converted_number


def stacked_bar_plot_single_sample(data, show_value, output):
    """
    Plot a single taxonomic data set.
    Require matplotlib and numpy
    :param data: The data set to plot
    :param show_value: Bool. True if you want to show values over each bar
    :return: Nothing but will open a plot window.
    """
    number_of_sample = len(data)
    taxon_name = []
    taxon_value = []
    for taxon in data:
        taxon_name.append(taxon)
        taxon_value.append(data[taxon])
    ind = np.arange(number_of_sample)
    width = 0.75    # could be optimized
    plot = ppl.bar(ind, taxon_value, width, color='g')
    ppl.ylabel("Proportion")
    ppl.title("Taxonomic repartition")
    locs, labels = mpl.xticks(ind+width/2., taxon_name)
    ppl.setp(labels, rotation=90)
    ppl.yticks(np.arange(0, max(taxon_value), max(taxon_value)/10))
    if mpl.__version__ >="1.3.1":
        ppl.tight_layout()
    if show_value:
        for pos, value in zip(ind, taxon_value):
            ppl.text(pos + 0.5*width, value, str(value)[0:6], ha='center', va='bottom', rotation=90)
    if output is None:
        ppl.show()
    else:
        ppl.savefig(output)


def stacked_bar_plot_multi_sample(data, show_value, output):
    if show_value:
        sys.stderr.write("Warning: Can't show values for stacked plot\n")
    taxon_values = {}
    current_index = 1
    sample_list = []
    for sample in data:
        sample_list.append(sample)
        for taxon in data[sample]:
            if taxon in taxon_values:
                taxon_values[taxon].append(data[sample][taxon])
            else:
                taxon_values[taxon] = []
                i = 1
                while i < current_index:
                    taxon_values[taxon].append(0)
                    i += 1
                taxon_values[taxon].append(data[sample][taxon])
        for taxon in taxon_values:
            while len(taxon_values[taxon]) < current_index:
                taxon_values[taxon].append(0)
        current_index += 1
    number_of_sample = len(sample_list)
    ind = np.arange(number_of_sample)
    width = 0.35
    previous_taxon = []
    taxon_names = []
    graph_list = []
    color_index = 0.0
    for taxon in taxon_values:
        taxon_names.append(taxon)


        if len(previous_taxon) == 0:
            graph_list.append(ppl.bar(ind, taxon_values[taxon], width, label=taxon,
                                      color=mpl.cm.Set1(1.*color_index/len(taxon_values))))
            previous_taxon = taxon_values[taxon]
            color_index += 1

        else:
            graph_list.append(ppl.bar(ind, taxon_values[taxon], width, bottom=previous_taxon, label=taxon,
                              color=mpl.cm.Paired(1.*color_index/len(taxon_values))))
            color_index += 1
            list_index = 0
            for j in taxon_values[taxon]:
                previous_taxon[list_index] += j
                list_index += 1
    fontP = FontProperties()
    fontP.set_size('small')
    ppl.ylabel("Proportion")
    ppl.title("Taxonomic repartition")
    ppl.yticks(np.arange(0, max(previous_taxon), 0.05))
    locs, labels = ppl.xticks(ind+width/2., sample_list)
    ppl.setp(labels, rotation=90)
    ppl.legend(loc='upper center', bbox_to_anchor=(0.5, -0.33), ncol=len(taxon_names)/6, fancybox=True, shadow=True, prop=fontP)
    if mpl.__version__ >="1.3.1":
        ppl.tight_layout()
    #if show_value:
    #    for pos, value in zip(ind, taxon_value):
    #        ppl.text(pos + 0.5*width, value, str(value)[0:6], ha='center', va='bottom', rotation=90)
    if output is None:
        ppl.show()
    else:
        ppl.savefig(output)

def compare_bar_plot_multi_sample(data, show_value, output):
    taxon_values = {}
    current_index = 1
    sample_list = []
    all_values = []
    taxon_names = []
    width = 0.02
    for sample in data:
        sample_list.append(sample)
        for taxon in data[sample]:
            if taxon in taxon_values:
                taxon_values[taxon].append(data[sample][taxon])
            else:
                taxon_values[taxon] = []
                i = 1
                while i < current_index:
                    taxon_values[taxon].append(0)
                    i += 1
                taxon_values[taxon].append(data[sample][taxon])
        for taxon in taxon_values:
            while len(taxon_values[taxon]) < current_index:
                taxon_values[taxon].append(0)
        current_index += 1
    number_of_sample = len(sample_list)
    ind = np.arange(number_of_sample)
    color_index = 0
    for taxon in taxon_values:
        taxon_names.append(taxon)
        all_values += taxon_values[taxon]
        ppl.bar(ind+(width*color_index), taxon_values[taxon], width, label=taxon,
                color=mpl.cm.Set1(1.*color_index/len(taxon_values)))
        color_index += 1

    fontP = FontProperties()
    fontP.set_size('small')
    ppl.ylabel("Proportion")
    ppl.title("Taxonomic repartition")
    ppl.yticks(np.arange(0, max(all_values), max(all_values)/10))
    locs, labels = ppl.xticks(ind+width/2., sample_list)
    ppl.setp(labels, rotation=90)
    ppl.legend(loc='upper center', ncol=len(taxon_names) % 6, bbox_to_anchor=(0.5, -0.33), fancybox=True, shadow=True,
               prop=fontP)
    if mpl.__version__ >="1.3.1":
        ppl.tight_layout()
    #if show_value:
    #    for pos, value in zip(ind, taxon_value):
    #        ppl.text(pos + 0.5*width, value, str(value)[0:6], ha='center', va='bottom', rotation=90)
    if output is None:
        ppl.show()
    else:
        ppl.savefig(output)






def read_taxonomy_multi_file(file_name, minimum_proportion):
    data = {}
    for taxon_file in open(file_name, 'r'):
        taxon_file_name = taxon_file.rstrip('\n')
        data[taxon_file_name] = read_taxonomy_file(taxon_file_name, minimum_proportion)
    return data


def read_db_multi_file(file_name, minimum_proportion):
    data = {}
    for taxon_file in open(file_name,'r'):
        taxon_file_name = taxon_file.rstrip('\n')
        data[taxon_file_name] = read_db_file(taxon_file_name, minimum_proportion)
    return data


#for axis rotation:
#http://stackoverflow.com/questions/10998621/rotate-axis-text-in-python-matplotlib
#for legends:
#http://matplotlib.org/users/legend_guide.html

if __name__ == "__main__":
    """
    New tool to plot Biological abundance directly from Ray output.
    Objective is to obtain a graph in less than a second
    Futur objective, a web page with possibly some parameters for sample exploration
    """
    parser = OptionParser(sys.argv[1:])
    args = parser.getArguments()
    if args["o"]is not None:
        mpl.use("Agg")
    global ppl
    import matplotlib.pyplot as ppl
    if sys.argv[1] != "plot":
        raise RuntimeError("The only option available in this script is 'plot'")
    data = {}
    if args["f"] is not None and args["multi"]==False:
        if args["t"] == "taxonomy":
            data = read_taxonomy_file(args["f"], args["m"])
        elif args["t"] == "db":
            data = read_db_file(args["f"])
        else:
            raise ValueError("The file type specified is invalid. Must be 'taxonomy' or 'db'")
        stacked_bar_plot_single_sample(data, args["value"], args["o"])
    if args["f"] is not None and args["multi"] is True:
        if args["t"] == "taxonomy":
            data = read_taxonomy_multi_file(args["f"], args["m"])
        if args["t"] == "db":
            data = read_db_multi_file(args["f"], args["m"])
        if args["stacked"] is True:
            stacked_bar_plot_multi_sample(data, args["value"], args["o"])
        else:
            compare_bar_plot_multi_sample(data, args["value"], args["o"])











