import sys,os
import numpy as np
import argparse

class Params:
	def __init__(self,folder = "data"):
		filename = "%sparam.00000.p00000"%folder
		self.d = {}
		with open(filename) as f:
			for line in f:
				(key, val) = line.split()
				self.d[key] = val				
	def get(self):
		return self.d
		
		
class ParamsInfo:
	def __init__(self,folder = "data"):
		filename = "%sparam.info"%folder
		self.d = {}
		with open(filename) as f:
			for line in f:
				(key, val) = line.split()				
				self.d[key] = float(val)
	def get(self):
		return self.d
		
class ParamsRun:
	def __init__(self,folder = "data"):
		filename = "%sparam.run"%folder
		self.d = {}
		with open(filename) as f:
			for line in f:				
				if line[0]!="#":
					(key, val) = line.split()					
					self.d[key] = val
	def get(self):
		return self.d
		
class StepAvg:
	"""
	step, aexp, z, max_level,max_rho,mean_xion,mean_T,max_T,stars,sfr
	"""
	def __init__(self,folder = "data"):
		filename = "%sstep.avg"%folder
	
		
		self.data= np.loadtxt(filename,skiprows=1,unpack=True)
	def get(self):
		return self.data
		
		
def printparam(folder):
	d = Params(folder).get()
	for key in d:
		print  key, d[key]

def comparParam(f1,f2):
	P1 = Params(f1).get()
	P2 = Params(f2).get()
	for key in P1:
		if not P1[key] == P2[key]:
			print key, P1[key], P2[key]

def splitPath(filepath):
	all_path = filepath.split("/") 
	filename = all_path[-1]
	folder = "/".join(all_path[:-1]) + "/"	

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


def readParamInfo(folder):
		return ParamsInfo(folder).get()

def readParamRun(folder):
		return ParamsRun(folder).get()

