import numpy as np

class Info:
	def __init__(self, folder):

		filename = "%s%s"%(folder,"data/param.info")
		f = open(filename)

		for line in f:
			if line[0]!="#":
				(key, val) = line.split()
				try :
					val=np.float64(val)
				except ValueError:
					pass
				setattr(self,key,val)

class Run:
	def __init__(self,folder):

		filename = "%s%s"%(folder,"SRC/param.run")
		f = open(filename)

		for line in f:
			if line[0]!="#":
				(key, val) = line.split()
				try :
					val=np.float64(val)
				except ValueError:
					pass
				setattr(self,key,val)

class Avg:
	def __init__(self,folder):

		filename = "%s%s"%(folder,"data/param.avg")
		f = open(filename)
		header=f.readline().split("\t")
		data= np.loadtxt(filename,skiprows=1,unpack=True)

		i=0
		for field in header:
			if (field !='')&(field !='\n') :				
				try:
					setattr(self,field,data[i])
				except IndexError:
					print(field)
					#print ("WARNING : Problem while reading param.avg")
					pass
				i+=1
				

class Field:
	def __init__(self,folder,field):
		filename = "%s%s%s"%(folder,"data/avg/",field)
		f = open(filename)
				
		data= np.loadtxt(filename,unpack=True)
		setattr(self,"a",data[0])
		setattr(self,"mean",data[1])
		setattr(self,"sigma",data[2])
		setattr(self,"min",data[3])
		setattr(self,"max",data[4])		


class Param:
	def __init__(self,folder):
		print("Reading param")
		self.run=Run(folder)
		self.avg=Avg(folder)
		self.info=Info(folder)
		
		
