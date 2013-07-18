#!/usr/bin/python


import sys;

def warning(statement):
	sys.stderr.write("Warning:"+str(statement)+"\n");

def error(statement):
	sys.stderr.write("Error:"+str(statement)+"\n");
	sys.exit(1);

