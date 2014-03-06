#!/usr/bin/python

"""
This script determines the X best matches from the DataBase used with RayMeta for each Contigs of an assembly.
Usage:
python FindContigsIdWithBiologicalAbundance.py BiologicalAbundanceDirectory numberOfBestMatch FileContainingThePathToContigsIdentification.tsvFiles(optional)
If a file with the directory of ContigsIdentification.tsv files  is given, only those files will be considered.
The program display only entries for wich there is a PL-value higher than 0. It can be disabled by modifying a fct in Biological....py
"""

import sys, os, subprocess
from Modules.BiologicalAbundanceContigObject import *
from Modules.Error import *


def testValidity():
	if(len(sys.argv)==3):
		try:
			directory=sys.argv[1];
			numberOfBestMatch=int(sys.argv[2]);
		except:
			print(__doc__);
			sys.exit(1);
	elif(len(sys.argv)==4):
        	try:
	                directory=sys.argv[1];
			pathToContigsIdFile=sys.argv[3];
			tempFile=open(pathToContigsIdFile, 'r');
	                tempFile.close();
		        numberOfBestMatch=int(sys.argv[2]);
		except:
			print(__doc__);
	                sys.exit(1);
	else:
		print(__doc__);
		sys.exit(1);
        try:   
                tempFile=open(directory+"/_DeNovoAssembly/Contigs.tsv", 'r');
                tempFile.close();
        except:
                error("Unable to open: "+directory+"/_DeNovoAssembly/Contigs.tsv");
                sys.exit(1);

def readPathsFile(pathToContigsIdFile):
	#TODO: verifier que si les fichiers ont ou non le ContigIdentifications.tsv a la fin.
	contigIdentificationsFiles=[];
	try:
		tempfile=open(pathToContigsIdFile,'r');
		tempfile.close();
	except:
		error("Unable to open: "+pathToContigsIdFile);
		sys.exit(1);
	for lines in open(pathToContigsIdFile, 'r'):
                contigIdentificationsFiles.append(lines[0:len(lines)-1]);
	return contigIdentificationsFiles;

def getPathsFromDirectory(directory):
	command="find "+directory+"| grep ContigIdentifications.tsv"
	files=subprocess.check_output(command, shell=True);
	contigIdentificationsFiles=files.split('\n');
	for entry in contigIdentificationsFiles:
		if(entry=="" or '\n' in entry):
			contigIdentificationsFiles.remove(entry);
	return contigIdentificationsFiles;

def readContigsTSVfile(pathToFile):
	directory=pathToFile+"/_DeNovoAssembly/Contigs.tsv";
        biologicalAbundanceContigs={};
        for lines in open(directory, 'r'):
                if(lines==""):
                        break;
                elif(lines[0]=="#"):
                        continue;
                else:
                        name=lines.split()[0];
                        length=int(lines.split()[2]);
                        coloredKmer=int(lines.split()[3]);
                        newBiologicalAbundanceContig=BiologicalAbundanceContigObject(name,length,coloredKmer);
                        #print(name);
			biologicalAbundanceContigs[name]=newBiologicalAbundanceContig;
	return biologicalAbundanceContigs;

def readContigIdentificationFiles(files, biologicalAbundanceContigs):
	for line in open(files, 'r'):
		if(line==""):
			break;
		if(line[0]=='#'):
			continue;
		else:
			try:
				#print(line);
				contigName=line.split()[0];
				#print("1");
				sequenceName=line.split('\t')[6];
				#print("2");
				sequenceLength=int(line.split('\t')[7]);
				#print("3");
				matches=int(line.split('\t')[8]);
				#print("4");
				biologicalAbundanceContigs[contigName].addNewContigIdentification(sequenceName, sequenceLength, matches);
				#print("5");
				#print(contigName);	
				#newBiologicalAbundanceContig.showTSV();
			except:
				message="File: "+files+" can not be read or do not correspond to the usual patern. It has been discarded";
				warning(message);

def showTSV(biologicalAbundanceContigs):
	header="Contig name"+'\t'+"Contig length in k-mer"+'\t'+"Contig colored k-mer"+'\t'+"Sequence name"+'\t'+"Sequence length in k-mer"+"\t"+"Matches in sequence"+'\t'+"pl-value";
	print(header);
	for contigs in biologicalAbundanceContigs:
		biologicalAbundanceContigs[contigs].showTSV(); 

if __name__=="__main__":
	testValidity();
	directory=sys.argv[1];
        numberOfBestMatch=int(sys.argv[2]);
	if(len(sys.argv)==3):
		pathToContigsIdFile="";
	else:
		pathToContigsIdFile=sys.argv[3];	
	#First step: read the Path file.
	contigIdentificationsFiles=[];
	if(pathToContigsIdFile!=""):
		contigIdentificationsFiles=readPathsFile(pathToContigsIdFile);
	else:
		contigIdentificationsFiles=getPathsFromDirectory(directory);
	#Second step: read the Contigs.tsv file from _DeNovoAssembly directory.
	biologicalAbundanceContigs={};
	biologicalAbundanceContigs=readContigsTSVfile(directory);
	#third step: For each ContigsIdentification.tsv file, read each line and put it at the right place.
	for files in contigIdentificationsFiles:
		try:
			readContigIdentificationFiles(files, biologicalAbundanceContigs);
		except:
			warning("Unable to read: "+ files);
	#Fourth step: Calculate PL-values.
	#print("contig-693000039 before PL values");
	#print(len(biologicalAbundanceContigs["contig-693000039"].contigIdentifications));
	for contigs in biologicalAbundanceContigs:
		biologicalAbundanceContigs[contigs].calculatePLvalues();
	#print("contig-693000039 after PL values");
        #print(len(biologicalAbundanceContigs["contig-693000039"].contigIdentifications));
	#new step: remove null PL-values entry
	#for contigs in biologicalAbundanceContigs:
	#	print(biologicalAbundanceContigs[contigs].contigIdentifications);
	#	biologicalAbundanceContigs[contigs].contigIdentifications=biologicalAbundanceContigs[contigs].removeNulIdentification();
	#	print(biologicalAbundanceContigs[contigs].contigIdentifications);
	#Fifth step: Select only the top ?? identification for each contigs
	for contigs in biologicalAbundanceContigs:
		if(len(biologicalAbundanceContigs[contigs].contigIdentifications)<=0):
			continue;
		else:
			biologicalAbundanceContigs[contigs].selectBestIdentifications(numberOfBestMatch);
	#biologicalAbundanceContigs["contig-693000039"].selectBestIdentifications(numberOfBestMatch);
	#print("contig-693000039 after top selection");
        #print(len(biologicalAbundanceContigs["contig-693000039"].contigIdentifications));
	#Last step: Write to stdout!
	showTSV(biologicalAbundanceContigs);
