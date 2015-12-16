import os
import numpy as np

def getnproc(path):	
	try:
		nproc = len(os.listdir(path))
	except FileNotFoundError:
		print("ERROR : can't determine nproc in \n%s"%path)
		pass
	return nproc
	
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

			#print(xmin,xmax,ymin,ymax,zmin,zmax)			
			#print(bound[4],zmin,bound[4]<=zmin)
			#print(bound[5],zmax,bound[5]>=zmax)


			#print(bound[0],xmin,bound[0]<=xmin)
			#print(bound[1],xmax,bound[1]>=xmax)
			#print (bound[0]<=xmin or bound[1]>=xmax)
			#print (bound[4]<=zmin and bound[5]>=zmax)
			#print (bound[4]<=zmin and bound[5]>=zmax)
			
			#if ( bound[0]<=xmin or bound[1]>=xmax or bound[2]<=ymin or  bound[3]>=ymax or bound[4]<=zmin or bound[5]>=zmax ):
			
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
