import os

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
		self.folder=folder
		self.data_folder=folder+"data/"
		
		for folder in os.listdir(self.data_folder):
			try:				
				stepnum=int(folder)							
			except ValueError:
				continue
				
			key="step_%05d"%stepnum
			val= step.Step(stepnum, self.data_folder)
			setattr(self,key,val)
		
		self.param=param.Param(self.folder)
		
class Step:
	def __init__(self,number,folder):
		"""			
		Create a step object			
		"""
	
		self.number=number
		self.folder=folder
				
		self.part=Fields(number,folder,0)
		self.star=Fields(number,folder,1)
		self.grid=Fields(number,folder,2)
		
		self.hop=hop.Hop(number,folder)


class Fields:
	def __init__(self, number,folder, sets_type):
		"""			
			Create a set of fields object
		"""
		
		self.number=number		
		self.sets_type	= sets_type
				
		if sets_type==0:
			self.type="star_"
		elif sets_type==1:
			self.type="part_"
		elif sets_type==2:
			self.type="grid_"
			
		path = "%s%05d/"%(folder,number)
		for cur_folder in  os.listdir(path):
			if self.type in cur_folder:
				key=cur_folder[5:].replace(".","_")
				val= field.Field(folder,number,cur_folder)
				setattr(self,key,val)
	
class Field():
	def __init__(self,runpath,stepnum,field):
		self.runpath=runpath     
		self.stepnum=stepnum
		self.num="{:05d}".format(self.stepnum)
		self.field=field
		
		self.field_folder="%s%s/%s/"%(self.runpath,self.num,self.field)			
		cur_field = field[5:]				
		self.filename="%s%s.%s"%(self.field_folder,cur_field,self.num)
				
	def read1proc(self,filename):
		with open(filename, "rb") as file:            
			N = np.fromfile(file, dtype=np.int32  ,count=1)[0]
			tsim = np.fromfile(file, dtype=np.float32  ,count=1)[0]        
			bound= np.fromfile(file, dtype=np.float32  ,count=6)

			bx= bound[0]>self.xmax or bound[1]<self.xmin
			by= bound[2]>self.ymax or bound[3]<self.ymin
			bz= bound[4]>self.zmax or bound[5]<self.zmin
			
			if bx or by or bz :            
				return 0,tsim,bound,[]
			else:
				data = np.fromfile(file, dtype=np.float32)
				return N,tsim,bound,data        
				
	def read(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1):
		self.xmin=xmin
		self.xmax=xmax
		self.ymin=ymin
		self.ymax=ymax
		self.zmin=zmin
		self.zmax=zmax      

		self.N=0
		self.data=[]
		self.bound=[]

		for proc in range(getnproc(self.field_folder)):            
			N1proc,tsim,bound1proc,data1proc=self.read1proc(self.filename + ".p"+ str(proc).zfill(5))
			self.N+=N1proc
			self.tsim=tsim            
			self.bound=np.append(self.bound,bound1proc)
			self.data=np.append(self.data,data1proc)        

