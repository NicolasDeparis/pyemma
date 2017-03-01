# coding: utf-8

import os
import numpy as np
import h5py

import pyemma.fof as fof
# import hop
# import optical_depth

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

        self._hdf5 = True

        for folder in np.sort(os.listdir(self._data_folder)):
            try:
                stepnum=int(folder)
            except ValueError:
                continue

            key="step_%05d"%stepnum
            val= Step(stepnum, self._data_folder)
            setattr(self,key,val)

            self.step_list.append(val)
        self.param=Param(self.folder)

class Step:
    def __init__(self,number,folder):
        """
        Step object
        Contain all the data associated to an output time step
        """

        self.n=number # snap number

        self.part=Fields(number,folder,"part_") #particles
        self.star=Fields(number,folder,"star_") #stars
        self.grid=Fields(number,folder,"grid_") #AMR grid

#        self.optical_depth=optical_depth.OpticalDepth() #Optical depth

        self.a=self.grid._get_a() # scale factor
        self.z=1./self.a -1 # redshift


#comment these line in case of probleme with halo finder
#        self.hop=hop.Hop(number,folder)
        self.fof=fof.Fof(folder,number)

class Fields:
    def __init__(self, number,folder, sets_type):
        """
        Create a set of field object
        """

        self._number=number
        self._type=sets_type

        path = "%s%05d/"%(folder,number)
        for cur_folder in  np.sort(os.listdir(path)):

            if not "h5" in cur_folder:
                continue
            if  os.path.isdir("%s%s"%(path,cur_folder)):
                continue

            if self._type in cur_folder:
                key=(cur_folder[5:].replace(".","_").replace("[","").replace("]",""))[:-9]
                val= Field(folder,number,cur_folder[:-9])
                setattr(self,key,val)

    def _get_a(self,force=0):
        """
        get the scale factor.
        """
        for i in self.__dict__ :
            if self.__dict__[i].__class__.__name__ == "Field":
                return self.__dict__[i]._get_a()

class Field():
    """
    Field object (this is the main data object)
    """

    def __init__(self,runpath,stepnum,field):
        self._runpath=runpath
        self._stepnum=stepnum
        self._field=field
        self._field_folder="%s%05d/"%(runpath,stepnum)
        self._filename="%s%s_%05d.h5"%(self._field_folder,field,stepnum)
        self._isloadded=False

    def __getattr__(self, name):
        """
        automatic call to read when asking for data
        """
        if name == 'data':
            self.read()
            # self.read_old()
            return self.__getattribute__(name)
        else:
            raise AttributeError

    def read(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1, force=0):
        """
        The main reader function
        """
        if not self._isloadded or force :
            #print("Reading %s"%self._field)

            f = h5py.File(self._filename, "r")
            self.data=f['data'][:]
            f.close()

            self._isloadded=True
        else:
            print("%s allready loaded, use force=1 to reload"%self._field)

    def _get_a(self):
        """
        read the scale factor
        """
        f = h5py.File(self._filename, "r")
        return f.attrs['a']

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

    def _read_old(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1, force=0):

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

# # Parameter files

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

class RunParam:
    """
    Reader for param.run
    """
    def __init__(self,folder):
        filename = "%s%s"%(folder,"param.run")
        with open(filename) as f:
            for line in f:
                if line[0]!="#":
                    (key, val) = line.split()
                    try :
                        val=np.float64(val)
                    except ValueError:
                        pass
                    setattr(self,key,val)

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

class Avg:
    """
    Reader for param.avg.cpu
    """
    # TODO add cpu/gpu switch
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

class Param:
    def __init__(self,folder):
        self.run=RunParam(folder)
        self.avg=Avg(folder)
        self.info=Info(folder)
