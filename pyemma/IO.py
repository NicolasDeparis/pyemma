import sys,os
import matplotlib.pylab as plt
import numpy as np
import argparse


def splitPath(filepath):
	all_path = filepath.split("/") 	
	filename = all_path[-1]
	folder = "/".join(all_path[:-1]) + "/"	
	return folder, filename

def splitPath2(filepath):
	all_path = filepath.split("/") 
	print "allpath", all_path
	filename = all_path[-1]
	folder = "/".join(all_path[:-3]) + "/"
	print folder, filename
	return folder, filename
	
def getNproc(filepath):
	folder, filename = splitPath(filepath)
	files = os.listdir(folder)
	nProc =0
	for file in files:
		if filename in file :
			nProc += 1
	return nProc
	
def getargs() :
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('-fo',   action="store",      default="data/",     help = "witch folder to use", dest = "folder")
	parser.add_argument('-fi',   action="append",     default=[],            help = "snapshot file(s)", dest = "files")
	parser.add_argument('-np',   action="store",      default=-1, type=int,  help = "number of procs used to generate the snap. only usefull to force it", dest = "nproc")
	parser.add_argument('-l'    ,action="store",      default=0 , type=int,  help = "level param of oct2grid", dest = "level")
	parser.add_argument('-field',action="append",     default=["field.d"] ,  help = "field param of oct2grid", dest = "field")

	args = parser.parse_args()
	return args
