# coding: utf-8

import os
import numpy as np
import h5py

import fof
# import hop

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

class Step:
    def __init__(self,number,folder):
        """
        Step object
        """

        self.part=Fields(number,folder,"part_")
        self.star=Fields(number,folder,"star_")
        self.grid=Fields(number,folder,"grid_")

#comment these line in case of probleme with halo finder
#       self.hop=hop.Hop(number,folder)
        self.fof=fof.Fof(folder,number)

        self.a=self.grid._get_a()

class Fields:
    def __init__(self, number,folder, sets_type):
        """
        Create a set of field object
        """

        self._number=number
        self._type=sets_type

        path = "%s%05d/"%(folder,number)
        for cur_folder in  os.listdir(path):

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
            return self.__getattribute__(name)
        else:
            raise AttributeError

    def read(self, xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1, force=0):
        """
        The main reader function
        """
        if not self._isloadded or force :
            print("reading %s"%self._field)

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
