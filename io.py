import os
import numpy as np

import hop
import param

def getnproc(path):	
	try:
		nproc = len(os.listdir(path))
	except FileNotFoundError:
		print("ERROR : can't determine nproc in \n%s"%path)
		pass
	return nproc

class Run:
	def __init__(self,folder):
		self._folder=folder
		self._data_folder=folder+"data/"
		
		for folder in os.listdir(self._data_folder):
			try:				
				stepnum=int(folder)							
			except ValueError:
				continue
				
			key="step_%05d"%stepnum
			val= Step(stepnum, self._data_folder)
			setattr(self,key,val)
		
		self.param=param.Param(self._folder)
		
class Step:
	def __init__(self,number,folder):
		"""			
		Create a step object			
		"""
	
		self._number=number
		self._folder=folder
				
		self.part=Fields(number,folder,0)
		self.star=Fields(number,folder,1)
		self.grid=Fields(number,folder,2)
		
		self.hop=hop.Hop(number,folder)
		
class Fields:
	def __init__(self, number,folder, sets_type):
		"""			
		Create a set of fields object
		"""
		
		self._number=number		
		self._sets_type	= sets_type
				
		if sets_type==0:
			self._type="star_"
		elif sets_type==1:
			self._type="part_"
		elif sets_type==2:
			self._type="grid_"
			
		path = "%s%05d/"%(folder,number)
		for cur_folder in  os.listdir(path):
			if self._type in cur_folder:
				key=cur_folder[5:].replace(".","_")
				val= Field(folder,number,cur_folder)
				setattr(self,key,val)
	
class Field():
	"""
	Field object
	"""
	def __init__(self,runpath,stepnum,field):
		self._runpath=runpath     
		self._stepnum=stepnum		
		self._field=field
		
		self._field_folder="%s%05d/%s/"%(runpath,stepnum,field)
		cur_field = field[5:]				
		self._filename="%s%s.%05d"%(self._field_folder,cur_field,stepnum)
				
	def _read1proc(self,filename):
		with open(filename, "rb") as file:            
			N = np.fromfile(file, dtype=np.int32  ,count=1)[0]
			tsim = np.fromfile(file, dtype=np.float32  ,count=1)[0]        
			bound= np.fromfile(file, dtype=np.float32  ,count=6)

			bx= bound[0]>self._xmax or bound[1]<self._xmin
			by= bound[2]>self._ymax or bound[3]<self._ymin
			bz= bound[4]>self._zmax or bound[5]<self._zmin
			
			if bx or by or bz :            
				return 0,tsim,bound,[]
			else:
				data = np.fromfile(file, dtype=np.float32)
				return N,tsim,bound,data        
				
	def read(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1):
		self._xmin=xmin
		self._xmax=xmax
		self._ymin=ymin
		self._ymax=ymax
		self._zmin=zmin
		self._zmax=zmax      

		self._N=0
		self.data=[]
		self._bound=[]

		for proc in range(getnproc(self._field_folder)):
			N1proc,tsim,bound1proc,data1proc=self._read1proc(self._filename + ".p"+ str(proc).zfill(5))
			self._N+=N1proc
			self._tsim=tsim            
			self._bound=np.append(self._bound,bound1proc)
			self.data=np.append(self.data,data1proc)        

