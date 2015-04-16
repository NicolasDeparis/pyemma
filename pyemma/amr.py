import numpy as np
import os,sys
import convert
import matplotlib.pylab as plt
import multiprocessing as mp
import scipy.ndimage

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
	dump = np.fromfile(file, dtype=np.int32    ,count=3)
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
	
	def __init__(self,file,field):
		
		
		convert.oct2cell(file,field=field)
		fileName=file.replace("grid.","alloct.") + "." + field

		
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
							
	def getN(self):
		return self.n
	def getPos(self):
		return self.x,self.y,self.z
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
	def min(self):
		return np.min(self.map)
	def histogram(self):
		plt.hist(self.map,50,log=True)
		plt.show()
	def lmin(self):
		return np.min(self.level)
	def lmax(self):
		return np.max(self.level)
		
	def plot(self):		
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


			for i in range(mask[0].shape[0]):
				grid[y[i],x[i]] += map_tmp[i]
			
			if (l<lmax):
				img= scipy.ndimage.zoom(grid,np.power(2,(lmax-l)),order=0)
				grid_full+= img
			else:
				grid_full+=grid

		plt.figure()
		plt.imshow(np.log10(grid_full),interpolation="none", origin='lower')
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
		
