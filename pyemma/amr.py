import numpy as np
import os,sys
import matplotlib.pylab as plt
import multiprocessing as mp
import scipy.ndimage

import IO
import struct
import convert

class array :
	def __init__(self,fileName):
		file = open(fileName, "rb")
		self.x,self.y,self.z = np.fromfile(file, dtype=np.int32    ,count=3)
		print self.x,self.y,self.z
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

def cube2log(fileName):
	DATA = array(fileName)
	x,y,z= DATA.getSize()
	print x,y,z
	a=DATA.geta()
	data = DATA.getData()
	print "read ok"
	
	f = open(fileName + ".log10", "wb")
	f.write(struct.pack('i', x))
	f.write(struct.pack('i', y))
	f.write(struct.pack('i', z))
	f.write(struct.pack('f', a))
	
	data = np.log10(data)
#	data *= 255/ (np.max(data)-np.min(data))
	
	for a in data:
		f.write(struct.pack('f', a)) 
	
	f.close()
	
def test_cube2log(fileName):
	DATA = cube(fileName)
	DATA = cube(fileName + ".log10")
	
	
def geta(fileName):
	file = open(fileName, "rb")
	#dump = np.fromfile(file, dtype=np.int32    ,count=3)
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
	
	
####################################################################################################
####################################################################################################


class allOct :
	
	def __init__(self,file,field="field.d",  xmin=0, xmax=-1, ymin=0, ymax=-1, zmin=0, zmax=-1,force=0):
		
		convert.oct2cell(file,field=field ,  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,force=force)
		fileName=file.replace("grid.","alloct.") + "." + field

		self.field = field
		
		NREAL=5
		type=np.float32

		file = open(fileName, "rb")
		self.n = np.fromfile(file,dtype=np.int32,count=1)[0]
		self.a = np.fromfile(file,dtype=np.float32,count=1)[0]		
		data = np.fromfile(file,dtype=type,count=self.n*5)
		file.close()

		self.x=np.zeros(self.n,dtype=type)
		self.y=np.zeros(self.n,dtype=type)
		self.z=np.zeros(self.n,dtype=type)
		self.level=np.zeros(self.n,dtype=type)
		self.map=np.zeros(self.n,dtype=type)

		for i in range(self.n):
			self.x[i]=data[i*NREAL+0]
			self.y[i]=data[i*NREAL+1]
			self.z[i]=data[i*NREAL+2]
			self.level[i]=data[i*NREAL+3]
			self.map[i]=data[i*NREAL+4]
			
		
	#	self.locatemax2()
		self.getMtot()	
	def getN(self):
		return self.n
	def geta(self):
		return self.a
	def getx(self):
		return self.x
	def gety(self):
		return self.y
	def getz(self):
		return self.z
	def getPos(self,n):
		return self.x[n],self.y[n],self.z[n]
	def getLevel(self):
		return self.level
	def getMap(self):
		return self.map
	def get(self):
		return self.getSize(),self.getLevel(),self.getMap()
	def mean(self):
		return np.mean(self.map)
	def max(self):
		return np.max(self.map)
	def locatemax(self):
		return np.argmax(self.map)	
	def locatemax2(self):
		arg = np.argmax(self.map)	
		return self.x[arg], self.y[arg], self.z[arg]
		
	def min(self):
		return np.min(self.map)
	def histogram(self):
		plt.hist(self.map,50,log=True)
		plt.show()
	def lmin(self):
		return np.min(self.level)
	def lmax(self):
		return np.max(self.level)
	def getMtot(self):		
		print np.sum(self.map * np.power(0.5, 3*self.level))
		
	def plot(self, log=0):		
		lmin=np.min(self.level)	
		lmax=np.max(self.level)
		#lmax=lmin
		grid_full=np.zeros((2**lmax,2**lmax))

		for l in range(lmin,int(lmax+1)):
			dx=np.power(2.,l)
			grid=np.zeros((dx,dx))		
			mask=np.where(self.level==l)
			print l,dx, mask[0].shape[0]
			
			x=np.int64(dx*self.x[mask])
			y=np.int64(dx*self.y[mask])
			map_tmp=self.map[mask]

			type ="mean"
			
			if type=="mean":
				dx = np.power(2,l)							
				for i in range(mask[0].shape[0]):
					grid[y[i],x[i]] += map_tmp[i]/dx
				
				if (l<lmax):
					grid_full+= scipy.ndimage.zoom(grid,np.power(2,(lmax-l)),order=0)
				else:
					grid_full+= grid
					
			if type=="max":
						
				for i in range(mask[0].shape[0]):
					grid[y[i],x[i]] = np.fmax(grid[y[i],x[i]],map_tmp[i])					
				
				if (l<lmax):
					tmp = scipy.ndimage.zoom(grid,np.power(2,(lmax-l)),order=0)		
					grid_full= np.fmax(tmp, grid_full)
				else:
					grid_full= np.fmax(grid, grid_full)
					



		plt.figure()
		if log:
			grid_full=np.log10(grid_full)
		plt.imshow(grid_full,interpolation="none", origin='lower')
		
		plt.colorbar()
		plt.show()
		
		
		
		
		
		"""
		def worker(mapall,iproc,xtot,ytot,nproc):
			N=self.n/nproc
			
			for i in range(iproc*N,(iproc+1)*N):
				curlev=self.level[i]
				if curlev<=lmax:
					
					curdl= int(2**(lmax-curlev))
					curx = int(self.x[i]*xtot) 
					cury = int(self.y[i]*ytot) 
					
					for lx in range(curdl):
						for ly in range(curdl):
							mapall[curx+lx][cury+ly]+= self.map[i]
		
		lmax=np.max(self.level)	
		lmax=8
		xtot=2**lmax
		ytot=2**lmax					
		
		tmp = np.ctypeslib.as_ctypes(np.zeros((xtot,ytot)))

		map_multi=mp.Array(tmp._type_,tmp, lock=False)

		nproc=6
		jobs = []
		for i in range(nproc):
			p = mp.Process(target=worker, args=(map_multi,i,xtot,ytot,nproc))
			jobs.append(p)
			p.start()
		for j in jobs:
			j.join()
	
		plt.imshow(map_multi, interpolation='nearest')
		plt.show()
		"""


####################################################################################################
####################################################################################################

def bound(xmin, xmax,ymin,ymax,zmin,zmax):
	if xmin<0:
		xmin=0	
	if xmax>1:
		xmax=1	
	if ymin<0:
		ymin=0	
	if ymax>1:
		ymax=1		
	if zmin<0:
		zmin=0	
	if zmax>1:
		zmax=1			
	return xmin,xmax,ymin,ymax,zmin,zmax

def default_min_max(xmin,xmax,ymin,ymax,zmin,zmax,N):
	if xmin == None:
		xmin = 0
	if ymin == None:
		ymin = 0
	if zmin == None:
		zmin = 0
	if xmax == None:
		xmax = N
	if ymax == None :
		ymax = N
	if zmax == None :
		zmax = N
	return xmin,xmax,ymin,ymax,zmin,zmax

def getSlice(filename,level,force=0, nproc=None, field="density", xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, log=True):
	
	xmin,xmax,ymin,ymax,zmin,zmax= default_min_max(xmin,xmax,ymin,ymax,zmin,zmax, 2**level)
	xmin,xmax,ymin,ymax,zmin,zmax= bound(xmin, xmax,ymin,ymax,zmin,zmax)
	
	if nproc==None:
		nproc=IO.getNproc(filename)
		
	diag=0
	if diag:
		convert.oct2grid(filename,level, force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, proj=0)
		data = cube(filename.replace("grid.","cube"+ str(level)+ ".") + "." + field)
		data = data.diagonal(1,0,1)
	else:
		convert.oct2grid(filename,level, force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, proj=3)
		data = cube(filename.replace("grid.","slice"+ str(level)+ ".") + "." + field)
		
	print "Z=%.4f"%data.getZ()
	print "a=%.4f"%data.geta()
	data = data.getData()
		
	if log:
		data = np.log10(data)
	

	#print data.shape
	#print data

	data = np.max(data,axis=0)
	#data = np.mean(data,axis=0)
	return data
	
def getCube(filename,level,force=0, nproc=None, field="density", xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None, log=False):
	
	xmin,xmax,ymin,ymax,zmin,zmax= default_min_max(xmin,xmax,ymin,ymax,zmin,zmax, 2**level)	
	xmin,xmax,ymin,ymax,zmin,zmax= bound(xmin, xmax,ymin,ymax,zmin,zmax)
	
	if nproc==None:
		nproc=IO.getNproc(filename)
	
	convert.oct2grid(filename,level, force, nproc, field, xmin, xmax, ymin, ymax, zmin, zmax, proj=0)
	data = cube(filename.replace("grid.","cube"+ str(level)+ ".") + "." + field)
	
	
	print "Z=%.4f"%data.getZ()
	print "a=%.4f"%data.geta()
	#data = data.getData()
		
	if log:
		data = np.log10(data)
		
	return data
		
####################################################################################################
####################################################################################################


		
class field :
	def __init__(self,fileName):	
		file = open(fileName, "rb")
		self.n = np.fromfile(file,dtype=np.int32,count=1)[0]		
		self.x = np.fromfile(file,dtype=np.float32,count=-1)
		file.close()

	def get(self):
		return self.n, self.x
	
	
class pos :	
	def __init__(self, folder, step_number):

		format_number = "{:05d}".format(step_number)
	
		field_name = "x"
		n,self.x = field("%s%s/%s.%s"%(folder,format_number,field_name,format_number)).get()			
		
		field_name = "y"
		n,self.y = field("%s%s/%s.%s"%(folder,format_number,field_name,format_number)).get()
		
		field_name = "z"
		n,self.z = field("%s%s/%s.%s"%(folder,format_number,field_name,format_number)).get()
		
		field_name = "l"
		n,self.l = field("%s%s/%s.%s"%(folder,format_number,field_name,format_number)).get()
				
		field_name = "field.d"
		n,self.map = field("%s%s/%s.%s"%(folder,format_number,field_name,format_number)).get()		
		
		if n != len(self.map):
			print n, len(self.map)
		
		lmin=np.min(self.l)		
		dx=np.int32(np.power(2.,lmin))
		

		h,xe,ye =np.histogram2d(self.y,self.z,bins=dx, weights=self.map)
				
		plt.imshow(np.log10(h),interpolation="none", origin='lower')
		plt.colorbar()
		
