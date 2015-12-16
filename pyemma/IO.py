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
