# coding: utf-8

import os
import numpy as np
import h5py

import fof
# import hop
import optical_depth
import movie


"""
TEST
"""

class Run:
    """
    Run object
    create a Param object and n Step object.
    look for subfolder of data/ with name convertible into int
    """

    def __init__(self,folder, hdf5=True):


        self.folder=folder
        self._data_folder=folder+"data/"
        self.step_list=[]
        self.nstep=0

        self._hdf5 = hdf5

        for folder in np.sort(os.listdir(self._data_folder)):
            try:
                stepnum=int(folder)
            except ValueError:
                continue

            key="step_%05d"%stepnum
            val= Step(stepnum, self._data_folder, hdf5)
            setattr(self,key,val)

            self.step_list.append(val)
            self.nstep+=1

        try:
            self.param=Param(self.folder)
        except FileNotFoundError:
            print('no param')
            pass

        try:
            self.movie=movie.Movie("%sdata/movie/"%self.folder)
        except FileNotFoundError:
            print('no movie')
            pass


class Step:
    def __init__(self,number,folder, hdf5=True):
        """
        Step object
        Contain all the data associated to an output time step
        """

        self.n=number # snap number

        self.part=Fields(number,folder,"part_",hdf5) #particles
        self.star=Fields(number,folder,"star_",hdf5) #stars
        self.grid=Fields(number,folder,"grid_",hdf5) #AMR grid


        self.optical_depth=optical_depth.OpticalDepth() #Optical depth

        # scale factor
        self.a=self.grid._get_a()
        if self.a is None:
            self.a=self.part._get_a()

        self.t=self.grid._get_t()
        if self.t is None:
            self.t=self.part._get_t()

        if self.a is not None:
            self.z=1./self.a -1 # redshift
        else :
            self.z=None

#comment these line in case of probleme with halo finder
#       self.hop=hop.Hop(number,folder)
        self.fof=fof.Fof(folder,number,self)

class Fields:
    def __init__(self, number,folder, sets_type, hdf5=True):
        """
        Create a set of field object
        """

        self._number=number
        self._type=sets_type

        path = "%s%05d/"%(folder,number)
        for cur_folder in  np.sort(os.listdir(path)):


            if hdf5:
                if not "h5" in cur_folder:
                    continue
                if  os.path.isdir("%s%s"%(path,cur_folder)):
                    continue

            if self._type in cur_folder:
                if hdf5:
                    key=(cur_folder[5:].replace(".","_").replace("[","").replace("]",""))[:-9]
                    val= Field(folder,number,cur_folder[:-9], hdf5)
                else:

                    key=(cur_folder[5:].replace(".","_").replace("[","").replace("]",""))
                    val= Field(folder,number,cur_folder, hdf5)

                setattr(self,key,val)



    def _get_a(self,force=0):
        """
        get the scale factor.
        """
        for i in self.__dict__ :
            if self.__dict__[i].__class__.__name__ == "Field":
                return self.__dict__[i]._get_a()

    def _get_t(self,force=0):
        """
        get the scale factor.
        """
        for i in self.__dict__ :
            if self.__dict__[i].__class__.__name__ == "Field":
                return self.__dict__[i]._get_t()

class Field():
    """
    Field object (this is the main reader object)
    """

    def __init__(self,runpath,stepnum,field,hdf5=True):
        self._runpath=runpath
        self._stepnum=stepnum
        self._field=field
        self._field_folder="%s%05d/"%(runpath,stepnum)
        self._isloadded=False
        self.hdf5=hdf5

        if hdf5:
            self._filename="%s%s_%05d.h5"%(self._field_folder,field,stepnum)
        else:
            self._filename="%s%s"%(self._field_folder,field)

    def __getattr__(self, name):
        """
        automatic call to read when asking for data
        """
        if name == 'data':
            self.read()
            return self.__getattribute__(name)
        else:
            raise AttributeError

    def _get_a(self):
        """
        read the scale factor
        """

        if self.hdf5:
            f = h5py.File(self._filename, "r")
            return f.attrs['a']
        else:
            f="%s/%s.%05d.p00000"%(self._filename,self._field[5:],self._stepnum)
            with open(f, "rb") as file:
                N = np.fromfile(file, dtype=np.int32  ,count=1)[0]
                if N==0:
                     return 0
                return np.fromfile(file, dtype=np.float32  ,count=1)[0]

    def _get_t(self):
        """
        read the physical time
        """
        if self.hdf5:
            try:
                f = h5py.File(self._filename, "r")
                return f.attrs['t']
            except KeyError:
                return 0
        else:
            return None

    def read(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1, force=0):
        """
        The main reader function
        """

        if self.hdf5:

            if not self._isloadded or force :
                print("reading %s"%self._field)

                f = h5py.File(self._filename, "r")
                self.data=f['data'][:]
                f.close()

                self._isloadded=True
            else:
                print("%s allready loaded, use force=1 to reload"%self._field)

        else:
             self._read_bin(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax, force=force)

    def _read1proc(self,filename):
        with open(filename, "rb") as file:
            N = np.fromfile(file, dtype=np.int32  ,count=1)[0]
            if N==0:
                 return 0,0,[],[]
            tsim = np.fromfile(file, dtype=np.float32  ,count=1)[0]

# # !!! Temporary fix !!!
            if "grid" in filename:
                bound= np.fromfile(file, dtype=np.float32  ,count=6)
            elif "part" in filename:
                bound=np.array([0,1,0,1,0,1])
            elif "star" in filename:
                bound=np.array([0,1,0,1,0,1])
            else:
                print("problem in reading field")
# # !!! end of temporary fix!!!

            # bound=np.array([0,1,0,1,0,1])

            bx= bound[0]>self._xmax or bound[1]<self._xmin
            by= bound[2]>self._ymax or bound[3]<self._ymin
            bz= bound[4]>self._zmax or bound[5]<self._zmin

            if bx or by or bz :
                return 0,tsim,bound,[]
            else:
                # print("reading %s "%filename)
                data = np.fromfile(file, dtype=np.float32)
                return N,tsim,bound,data


    def _getnproc(self):

        self.nproc=0
        try:
            self.files = os.listdir(self._filename)
            for file in self.files:
                if self._field[5:] in file :
                    self.nproc += 1
        except FileNotFoundError:
            pass

    def _read_bin(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1, force=0):

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

            self._getnproc()

            for proc in range(self.nproc):
                cur_name = "%s/%s.%05d.p%s"%(self._filename, self._field[5:] , self._stepnum, str(proc).zfill(5))
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
        filename = "%s%s"%(folder,"data/param.run")
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
