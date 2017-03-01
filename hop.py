# coding: utf-8

import os
import numpy as np

def getnproc(path):
    try:
        nproc = len( os.listdir(path+os.listdir(path)[0]) )
    except FileNotFoundError:
        print("ERROR : can't determine nproc")
        nproc=0
    return nproc

class Hop:
    def __init__(self,stepnum,folder):

        self.exec_folder  = "utils/hop/"

        self.number=stepnum
#         self.folder="%s/%05d"%(folder,stepnum)
        self.folder=folder

        self.name="hop.%05d"%stepnum
        self.path=("%s%s"%(folder,self.name))

        self.nproc=getnproc("%s%05d/"%(folder,stepnum))

        self.den_name = self.path + ".den"
        self.tag_name = self.path + ".tag"
        self.den_isloaded=False
        self.tag_isloaded=False

########################################################################

    def genDen(self):
        commande="./" + self.exec_folder + "hop -in " + self.folder + " -o " + self.path +" -p 1 -nf " + str(self.nproc) + " -step " + str(self.number)
        print ("executing %s"%commande)
        os.system(commande)

    def readDen(self):
        print ("Reading file", self.den_name)
        with open(self.den_name, "rb") as file:
            self.N = np.fromfile(file, dtype=np.int32   ,count=1)[0]
            self.den = np.fromfile(file, dtype=np.float32 ,count=self.N)
            print ("Read OK")

    def getDen(self,force=0):
        if (not os.path.isfile(self.den_name) or force):
            self.genDen()
        self.readDen()
        self.den_isloaded=True

    def genTag(self):
        commande =  "./" + self.exec_folder + "regroup -root %s -douter 80. -dsaddle 200. -dpeak 240. -f77 -o %s"%(self.path,self.path)
        print ("executing %s"%commande)
        os.system(commande)

    def readTag(self):
        print ("Reading file", self.tag_name)
        with open(self.tag_name, "rb") as file:
            dummy = np.fromfile(file, dtype=np.int32   ,count=1)
            self.npart = np.fromfile(file, dtype=np.int32   ,count=1)[0]
            self.ngrp = np.fromfile(file, dtype=np.int32   ,count=1)[0]
            dummy = np.fromfile(file, dtype=np.int32   ,count=1)
            dummy = np.fromfile(file, dtype=np.int32   ,count=1)
            self.tag = np.fromfile(file, dtype=np.int32   ,count=self.npart)
            dummy = np.fromfile(file, dtype=np.int32   ,count=1)
            print ("Read OK")

    def readNGRP(self):
        file = open(self.tag_name, "rb")
        dummy 	= np.fromfile(file, dtype=np.int32   ,count=1)
        self.npart = np.fromfile(file, dtype=np.int32   ,count=1)[0]
        self.ngrp = np.fromfile(file, dtype=np.int32   ,count=1)[0]
        file.close()

    def getTag(self,force=0):
        if (not os.path.isfile(self.tag_name) or force):
            self.genTag()
        self.readTag()
        self.tag_isloaded=True

########################################################################

    def get(self,force=0):
        self.getDen(force)
        self.getTag(force)
