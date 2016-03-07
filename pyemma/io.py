
# coding: utf-8

# In[1]:

import os
import numpy as np
import h5py


# In[16]:

def getnproc(path):
    """
    count the number of files in path
    """
    try:        
        nproc = len(os.listdir(path))
    except FileNotFoundError:
        print("ERROR : can't determine nproc in \n%s"%path)    
    return nproc


# In[2]:

class Run:
    """
    Run object 
    create a Param object and n Step object.
    look for subfolder of data/ with name convertible into int
    """

    def __init__(self,folder):
        self.folder=folder
        self._data_folder=folder+"data/"
        self.step_list=[]
        
        for folder in os.listdir(self._data_folder):
            try:
                stepnum=int(folder)
            except ValueError:
                continue
                    
            key="step_%05d"%stepnum
            val= Step(stepnum, self._data_folder)
            setattr(self,key,val)

            self.step_list.append(val)
        self.param=Param(self.folder)


# In[4]:

class Step:
    def __init__(self,number,folder ):
        """
        Step object
        """

        self.number=number
        self.folder=folder
        self.a=0.15

        self.part=Fields_hdf5(number,folder,0)
        self.star=Fields_hdf5(number,folder,1)
        self.grid=Fields_hdf5(number,folder,2)
  
#       import hop      
#       self.hop=hop.Hop(number,folder)

        import fof     
        self.fof=fof.Fof(folder,number)


# In[4]:

class Fields_bkp:
    def __init__(self, number,folder, sets_type):
        """
        Create a set of field object
        """

        self._number=number
        self._sets_type = sets_type

        if sets_type==0:
            self._type="part_"
        elif sets_type==1:
            self._type="star_"
        elif sets_type==2:
            self._type="grid_"

        path = "%s%05d/"%(folder,number)
        for cur_folder in  os.listdir(path):
            if self._type in cur_folder:
                key=cur_folder[5:].replace(".","_")
                key=key.replace("[","")
                key=key.replace("]","")
                val= Field(folder,number,cur_folder)
                setattr(self,key,val)
    
    def read(self,force=0):
        """
        read all field in fields.        
        """        
        for i in self.__dict__ :            
            if self.__dict__[i].__class__.__name__ == "Field":                
                self.__dict__[i].read(force)


# In[5]:

class Field_bkp():
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
        self._isloadded=False
        
    def __getattr__(self, name):                            
        if name == 'data':
            self.read()
            return self.__getattribute__(name)
        else:        
            raise AttributeError

    def _read1proc(self,filename):
        with open(filename, "rb") as file:                                    
            N = np.fromfile(file, dtype=np.int32  ,count=1)[0]
            if N==0:
                 return 0,0,[],[]
            tsim = np.fromfile(file, dtype=np.float32  ,count=1)[0]
# !!! Temporary fix !!!
            if "grid" in filename:    
                bound= np.fromfile(file, dtype=np.float32  ,count=6)
            elif "part" in filename:
                bound=np.array([0,1,0,1,0,1])
            elif "star" in filename:        
                bound=np.array([0,1,0,1,0,1])
            else:
                print("problem in reading field")
# !!! end of temporary fix!!!

            bx= bound[0]>self._xmax or bound[1]<self._xmin
            by= bound[2]>self._ymax or bound[3]<self._ymin
            bz= bound[4]>self._zmax or bound[5]<self._zmin
            
            if bx or by or bz :            
                return 0,tsim,bound,[]
            else:
                data = np.fromfile(file, dtype=np.float32)                
                return N,tsim,bound,data

    def read(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1, force=0):
        
        if not self._isloadded or force :
            print("reading %s"%self._field)
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
                cur_name = self._filename + ".p"+ str(proc).zfill(5)
                N1proc,tsim,bound1proc,data1proc=self._read1proc(cur_name)
                self._N+=N1proc
                self._tsim=tsim
                self._bound=np.append(self._bound,bound1proc)
                self.data=np.append(self.data,data1proc)
                
            self._isloadded=True
        else:
            print("%s allready loaded, use force=1 to reload"%self._field)


# In[4]:

class Fields:
    def __init__(self, number,folder, sets_type):
        """
        Create a set of field object
        """

        self._number=number
        self._sets_type = sets_type

        if sets_type==0:
            self._type="part_"
        elif sets_type==1:
            self._type="star_"
        elif sets_type==2:
            self._type="grid_"

        path = "%s%05d/"%(folder,number)
        for cur_folder in  os.listdir(path):
            if  os.path.isdir("%s%s"%(path,cur_folder)):
                continue
            if self._type in cur_folder:
                key=cur_folder[5:].replace(".","_").replace("[","").replace("]","")
                key=key[:-6]                
                val= Field(folder,number,cur_folder[:-6])
                setattr(self,key,val)
    
    def read(self,force=0):
        """
        read all field in fields.        
        """        
        for i in self.__dict__ :            
            if self.__dict__[i].__class__.__name__ == "Field":                
                self.__dict__[i].read(force)


# In[5]:

class Field():
    """
    Field object
    """
    def __init__(self,runpath,stepnum,field):
        self._runpath=runpath     
        self._stepnum=stepnum
        self._field=field

        self._field_folder="%s%05d/"%(runpath,stepnum)
        cur_field = field[5:]
        self._filename="%s%s.%05d"%(self._field_folder,field,stepnum)
        self._isloadded=False
        
    def __getattr__(self, name):                            
        if name == 'data':
            self.read()
            return self.__getattribute__(name)
        else:        
            raise AttributeError

    def read(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1, force=0):

        if not self._isloadded or force :
            print("reading %s"%self._field)
            with open(self._filename, "rb") as file:
                self.N = np.fromfile(file, dtype=np.int32  ,count=1)[0]
                self.nproc = np.fromfile(file, dtype=np.int32  ,count=1)[0]
                self.tsim = np.fromfile(file, dtype=np.float32  ,count=1)[0]
                self.ncells=np.fromfile(file, dtype=np.int32  ,count=self.nproc)
                self.bounds=np.fromfile(file, dtype=np.float32  ,count=6*self.nproc)
                self.data = np.fromfile(file, dtype=np.float32)

            self._isloadded=True
        else:
            print("%s allready loaded, use force=1 to reload"%self._field)


# In[ ]:

class Fields_hdf5:
    def __init__(self, number,folder, sets_type):
        """
        Create a set of field object
        """

        self._number=number
        self._sets_type = sets_type

        if sets_type==0:
            self._type="part_"
        elif sets_type==1:
            self._type="star_"
        elif sets_type==2:
            self._type="grid_"

        path = "%s%05d/"%(folder,number)
        for cur_folder in  os.listdir(path):
            if  os.path.isdir("%s%s"%(path,cur_folder)):
                continue
            if self._type in cur_folder:
                key=cur_folder[5:].replace(".","_").replace("[","").replace("]","")
                key=key[:-9]                
                val= Field_hdf5(folder,number,cur_folder[:-9])
                setattr(self,key,val)
    
    def read(self,force=0):
        """
        read all field in fields.        
        """        
        for i in self.__dict__ :            
            if self.__dict__[i].__class__.__name__ == "Field":                
                self.__dict__[i].read(force)


# In[ ]:

class Field_hdf5():
    """
    Field object
    """
    def __init__(self,runpath,stepnum,field):
        self._runpath=runpath     
        self._stepnum=stepnum
        self._field=field

        self._field_folder="%s%05d/"%(runpath,stepnum)
        cur_field = field[5:]
        self._filename="%s%s_%05d.h5"%(self._field_folder,field,stepnum)
        self._isloadded=False
        
    def __getattr__(self, name):                            
        if name == 'data':
            self.read()
            return self.__getattribute__(name)
        else:        
            raise AttributeError

    def read(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1, force=0):

        if not self._isloadded or force :
            print("reading %s"%self._field)
            
            f = h5py.File(self._filename, "r")
            self.data=f['data'][:]
            
            self._isloadded=True
        else:
            print("%s allready loaded, use force=1 to reload"%self._field)


# # Parameter files

# In[12]:

class Info:
    """
    Reader for param.info
    """
    def __init__(self, folder):
        filename = "%s%s"%(folder,"data/param.info")
        with open(filename) as f:

            for line in f:
                if line[0]!="#":
                    (key, val) = line.split()
                    try :
                        val=np.float64(val)
                    except ValueError:
                        pass
                    setattr(self,key,val)


# In[11]:

class RunParam:
    """
    Reader for param.run
    """
    def __init__(self,folder):
        filename = "%s%s"%(folder,"SRC/param.run")
        with open(filename) as f:
            for line in f:
                if line[0]!="#":
                    (key, val) = line.split()
                    try :
                        val=np.float64(val)
                    except ValueError:
                        pass
                    setattr(self,key,val)


# In[5]:

class FieldAvg:
    """
    Reader for avg/field.avg
    """
    def __init__(self,folder,field):
        filename = "%s%s%s"%(folder,"data/avg/",field)
        with open(filename) as f:
            data= np.loadtxt(filename,unpack=True)
        setattr(self,"a",data[0])
        setattr(self,"mean",data[1])
        setattr(self,"sigma",data[2])
        setattr(self,"min",data[3])
        setattr(self,"max",data[4])


# In[1]:

class Avg:
    """
    Reader for param.avg
    """
    def __init__(self,folder):
        filename = "%s%s"%(folder,"data/param.avg.cpu")
        with open(filename) as f:
            header=f.readline().split("\t")
            data= np.loadtxt(filename,skiprows=1,unpack=True)

            i=0
            for field in header:
                if (field !='')&(field !='\n') :
                    try:
                        setattr(self,field,data[i])
                    except IndexError:
                        pass
                    i+=1

            for cur_folder in  os.listdir(folder+"data/avg/"):
                try:
                    key = cur_folder.replace(".","_")
                    val = FieldAvg(folder,cur_folder)
                    setattr(self,key,val)
                except IndexError:
                    pass


# In[7]:

class Param:
    def __init__(self,folder):
        self.run=RunParam(folder)
        self.avg=Avg(folder)
        self.info=Info(folder)

