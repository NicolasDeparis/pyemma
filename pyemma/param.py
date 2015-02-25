import sys,os
import numpy as np
import matplotlib.pylab as plt

"""
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
"""

class Param:
	def __init__(self,folder = "data"):
		self.info=ParamsInfo(folder)
		self.run=ParamsRun(folder)
		self.avg=ParamsAvg(folder)
		
class ParamsInfo:
	def __init__(self,folder = "data"):
		filename = "%sparam.info"%folder
		self.d = {}
		with open(filename) as f:
			for line in f:				
				if line[0]!="#":		
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
		
class ParamsAvg:
	"""
	step, aexp, z, max_level,max_rho,mean_xion,mean_T,max_T,stars,sfr
	"""
	def __init__(self,folder = "data"):
		self.folder=folder
		filename = "%sparam.avg"%folder
		data= np.loadtxt(filename,skiprows=1,unpack=True)
		
		self.data = {}
		self.data["step"] = data[0]
		self.data["aexp"] = data[1]
		self.data["z"] = data[2]
		self.data["max_level"] = data[3]
		self.data["max_rho"] = data[4]
		self.data["mean_xion"] = data[5]
		self.data["mean_T"] = data[6]
		self.data["max_T"] = data[7]
		self.data["stars"] = data[8]

		
	def get(self):
		return self.data
	def evol(self,x_field,y_field, log=0,**kwargs):
		x=self.data[x_field]
		y=self.data[y_field]
		if log :
			plt.loglog(x,y,**kwargs)
		else:
			plt.plot(x,y, **kwargs)	
		plt.legend()
		plt.show()
		
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


