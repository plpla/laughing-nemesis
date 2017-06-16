__author__='pier-luc'

import glob
from datetime import date
import sys
allocation = "nne-790-ae"
script_path = "/rap/nne-790-ab/projects/pplante/Git_clone/laughing-nemesis/FindLastCommonAncester.py"
taxonomy_file = "/rap/nne-790-ab/projects/pplante/Taxonomy/TreeOfLife-Edges.tsv"
outpath = "/rap/nne-790-ab/projects/pplante/FromageSteveLabrie/Nemesis"

def create_new_sub_file(file_out, job_name):
    fo = open(file_out, 'w')
    ddate = str(date.today())
    fo.write("#!/bin/bash\n#PBS -N "+job_name+"\n#PBS -o "+job_name+".stdout\n" +
            "#PBS -e "+job_name+".stderr\n#PBS -A "+allocation+"\n#PBS -l walltime=8:00:00\n" +
            '#PBS -l nodes=1:ppn=8\n#PBS -q default\ncd "${PBS_O_WORKDIR}"\n\n')
    #fo.write("module load apps/python/2.7.5\n\n")
    fo.close()

def prepare_submissions_files(directory, prefix, job_name):
    sample_list = glob.glob(directory+"/*/BiologicalAbundances/")
    #print(sample_list)
    postfix = 1
    index = 0
    print("You can now submit theses files on Colosse:")
    for sample in sample_list:
        if index % 8 == 0:
            file_name = prefix+"_"+str(postfix)+".sh"
            print(file_name)
            create_new_sub_file(file_name, job_name+"_"+str(postfix))
            postfix += 1
        fo = open(file_name, 'a')
        fo.write("python "+script_path+" lca -d "+sample+" -t "+taxonomy_file+" >"+outpath+sample.split('/')[-3]+"_nemesis.tsv \n\n")
        #print("python "+script_path+" -d "+sample+" -t "+taxonomy_file+">"+sample+"_nemesis.tsv &\n")
        index += 1
        fo.close()
    fo = open(file_name, 'a')
    fo.write("\nwait\n\n")
    fo.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("To use: python Prepare_for_moab.py SAMPLE_DIRECTORY OUTPUT_PREFIX\n where the sample directory contains directories from RayMeta assemblies.\n")
        sys.exit(0)
    directory = sys.argv[1]
    prefix = sys.argv[2]
    job_name = "Nemesis"+str(date.today())
    prepare_submissions_files(directory, prefix, job_name)






