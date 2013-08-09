#!/usr/bin/python

import sys;


def isValid(file):
	valid=1;
	try:
        	a=open(file, 'r');
		a.close();
	except:
		valid=0;
	return valid;

def coutLines(file):
	lines=0;
	try:
		for i in open(file):
			lines+=1;
	except:
		lines=0;
	if lines==0:
		Error.error("File "+file+" is empty");
	return lines;
