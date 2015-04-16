import sys,os
import numpy as np
import matplotlib.pylab as plt

class Param:
	def __init__(self,folder = "data",stepnum=10):
		self.info=ParamsInfo(stepnum,folder)
		self.run=ParamsRun(stepnum,folder)
		self.avg=ParamsAvg(stepnum,folder)
		
class ParamsInfo:
	def __init__(self, stepnum,folder = "data"):
		paramfile="param.info"
		try :
			filename = "%s%s%s%s"%(folder,"{:05d}/".format(int(stepnum)),"grid/",paramfile)
			f = open(filename)		
		except  IOError:
			filename = "%s%s%s%s"%(folder,"{:05d}/".format(int(stepnum)),"part/",paramfile)
			f = open(filename)

		self.d = {}
		for line in f:				
			if line[0]!="#":		
				(key, val) = line.split()				
				self.d[key] = float(val)
	def get(self):
		return self.d
		
class ParamsRun:
	def __init__(self, stepnum,folder = "data"):
		paramfile="param.run"
		try :
			filename = "%s%s%s%s"%(folder,"{:05d}/".format(int(stepnum)),"grid/",paramfile)
			f = open(filename)		
		except  IOError:
			filename = "%s%s%s%s"%(folder,"{:05d}/".format(int(stepnum)),"part/",paramfile)
			f = open(filename)

		self.d = {}
		for line in f:				
			if line[0]!="#":
				(key, val) = line.split()					
				self.d[key] = val
	def get(self):
		return self.d
		
class ParamsAvg:
	"""
	step, aexp, z, dt, max_level,max_rho,mean_xion,mean_T,max_T,stars,sn
	"""
	def __init__(self, stepnum=0,folder = "data"):
		self.folder=folder
		
		paramfile="param.avg"
		try :
			filename = "%s%s%s%s"%(folder,"{:05d}/".format(int(stepnum)),"grid/",paramfile)
			f = open(filename)		
		except  IOError:
			filename = "%s%s%s%s"%(folder,"{:05d}/".format(int(stepnum)),"part/",paramfile)
			f = open(filename)

		data= np.loadtxt(filename,skiprows=1,unpack=True)
	
		self.data = {}
		self.data["step"] = data[0]
		self.data["aexp"] = data[1]
		self.data["z"] = data[2]
		self.data["dt"] = data[3]
		self.data["max_level"] = data[4]
		self.data["max_rho"] = data[5]
		self.data["mean_xion"] = data[6]
		self.data["mean_T"] = data[6]
		self.data["max_T"] = data[8]
		self.data["stars"] = data[9]
		self.data["SN"] = data[10]
	#	self.data["src"] = data[11]
		
	def get(self):
		return self.data
		
	def getxy(self,x_field,y_field,**kwargs):
		x=self.data[x_field]
		y=self.data[y_field]
		return x,y
		
	def evol(self,x_field,y_field, log=0,**kwargs):
		x=self.data[x_field]
		y=self.data[y_field]
		if log :
			plt.loglog(x,y,**kwargs)
		else:
			plt.plot(x,y, **kwargs)	
			
		plt.xlabel(x_field)
		plt.ylabel(y_field)
		plt.legend()
		
"""		
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
"""

def readParamInfo(folder):
		return ParamsInfo(folder).get()

def readParamRun(folder):
		return ParamsRun(folder).get()


