import sys,os
import numpy as np
import matplotlib.pylab as plt

import physique 
import sfr

class Param:
	def __init__(self,folder = "data",stepnum=10):
		self.info=ParamsInfo(stepnum,folder)
		self.run=ParamsRun(stepnum,folder)
	#	self.avg=ParamsAvg(stepnum,folder)
		
class ParamsInfo:
	def __init__(self, stepnum,folder = "data"):
		paramfile="param.info"
		try :
			filename = "%s%s"%(folder,paramfile)
			f = open(filename)		
		except  IOError:
			filename = "%s%s%s"%(folder,"../",paramfile)
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
			filename = "%s%s%s"%(folder,"../",paramfile)
			f = open(filename)		
		except  IOError:
			filename = "%s%s%s"%(folder,"../",paramfile)
			f = open(filename)

		self.d = {}
		for line in f:		
			
			if line[0]!="#":
				try :
					(key, val) = line.split()					
					self.d[key] = val
				except  IOError:		
					print line		
				
	def get(self):
		return self.d
		
class ParamsAvg:
	"""
	step, aexp, z, dt, max_level,max_rho,mean_xion,mean_T,max_T,stars,sn
	step	aexp		z		t_[yrs]		dt		max_level	max_rho		mean_xion	mean_T		max_T		stars		SN		src
	"""
	def __init__(self, stepnum=0,folder = "data"):
		self.folder=folder
		self.stepnum=stepnum
		paramfile="param.avg"
		try :
			filename = "%s%s"%(folder,paramfile)
			f = open(filename)		
		except  IOError:
			filename = "%s%s%s"%(folder,"../",paramfile)
			f = open(filename)

		data= np.loadtxt(filename,skiprows=1,unpack=True)
	
		self.data = {}
		self.data["step"] = data[0]
		self.data["aexp"] = data[1]
		self.data["z"] = data[2]
		self.data["t_[yrs]"] = data[3]
		self.data["dt"] = data[4]
		self.data["max_level"] = data[5]
		self.data["max_rho"] = data[6]
		self.data["mean_xion"] = data[7]
		self.data["mean_T"] = data[8]
		self.data["max_T"] = data[9]
		self.data["Nstars"] = data[10]
		self.data["Mstars"] = data[11]
		self.data["SN"] = data[12]
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
	
	def plot(self):
		f, axarr = plt.subplots(3,2, sharex=True)
		x=self.data["z"]
		
		y=self.data["max_level"]
		axarr[0,0].plot(x, y)
		axarr[0,0].set_title('level max')
		
		y=self.data["max_rho"]
		axarr[1,0].semilogy(x, y)
		axarr[1,0].set_title('rho max')
		
		y=self.data["mean_xion"]
		axarr[2,0].semilogy(x, y, label='xion')
		axarr[2,0].semilogy(x, 1.-y, label='xneu')
		axarr[2,0].set_title('ionisation')
		axarr[2,0].legend()
		
		y=self.data["mean_T"]
		axarr[0,1].semilogy(x, y)
		axarr[0,1].set_title('temp mean')
		
		y=self.data["max_T"]
		axarr[1,1].semilogy(x, y)
		axarr[1,1].set_title('temp max')
		
		y=self.data["stars"]
		p=ParamsInfo(self.stepnum,self.folder).get()
		y=np.diff(y)
		y*=p["mass_res_star"]
		y/=(p["box_size_Mpc/h"]/p["H0"]*100) **3
		y/=np.diff(physique.a2t(self.data["aexp"]))
		sfr.observation()
		print "instant sfr", y[-1:]
		
		x = 	(x[1:]+x[:-1])/2		

		#axarr[2,1].semilogy(x, y)
		axarr[2,1].set_title('stars')
		

				
		plt.xlim(5,20)
		
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


