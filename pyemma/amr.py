import numpy as np
import os
import convert

class array :
	def __init__(self,fileName):
		file = open(fileName, "rb")

		self.x,self.y,self.z = np.fromfile(file, dtype=np.int32    ,count=3)
		self.a = np.fromfile(file, dtype=np.float32  ,count=1)[0]
		self.data = np.fromfile(file, dtype=np.float32  ,count=self.x*self.y*self.z)
		file.close()
		
	def getData(self):
		return self.data
	def getN(self):
		return self.x*self.y*self.z
	def getSize(self):
		return self.z,self.y,self.x
	def geta(self):
		return self.a
	def getZ(self):
		return 1.0/self.a-1.0

class cube :
	def __init__(self,fileName):

		self.array = array(fileName)			
		self.data = np.reshape(self.array.getData(),  (self.array.getSize()) ) 
		
	def getData(self):
		return self.data
	def geta(self):
		return self.array.geta()
	def getSize(self):
		return self.array.getSize()
	def getZ(self):
		return 1.0/self.array.geta()-1.0
			
def geta(fileName):
	file = open(fileName, "rb")
	a = np.fromfile(file, dtype=np.float64  ,count=1)[0]
	file.close()
	return a
	
def printz(fileName):
	file = open(fileName, "rb")
	dump = np.fromfile(file, dtype=np.int32    ,count=3)
	a = np.fromfile(file, dtype=np.float32  ,count=1)[0]
	print "z = %s"%(1./a-1.)

	
def cubename(filename,level,field):
	return filename.replace("grid","cube"+ str(level)) + "." + str(field)

def gettmap(folder = "data/"):
	tmap = []
	for file in np.sort(os.listdir(folder)):
		if "grid." in file and ".p00000" in file: 
			file= "%s%s"%(folder,file)
			t = geta(file)
			tmap.append(t)
	return np.array(tmap)


def getXion(filename,level,force=0):
	
	convert.oct2grid(filename,level,force=force, field="rfield.xion")
	
	DATA = cube(cubename(filename,level,"rfield.xion"))
	data = DATA.getData()
	
	print DATA.getZ()
	print np.mean(data)
