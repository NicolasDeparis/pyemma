import numpy as np
import sys,os

import IO

def genFiles(step):

	filename = step.part_path

	folder  = "utils/hop/"

	outname = filename[:-10]  + "halo" + filename[-6:]

	genDen(folder,step)
	genTag(folder,step)
	#genPos(folder, outname)
	
########################################################################
def genDen(folder,step):
	
	filename = step.part_path
	outname = step.halo_path

	try :
		open(outname + ".den", "rb")
	except  IOError:
		nproc = step.nproc
		commande="./" + folder + "hop -in " + filename + " -o " + outname +" -p 1 -nf " + str(nproc)
		print commande
		os.system(commande)
		
def readDen(filename, step):

	denname = filename[:-10]  + "halo" + filename[-6:] + ".den"
	
	try :
		open(denname,"rb")
	except  IOError:
		genFiles(step)
		
	print "Reading file", denname
	file = open(denname, "rb")
	N = np.fromfile(file, dtype=np.int32   ,count=1)[0]
	den = np.fromfile(file, dtype=np.float32 ,count=N)
	file.close()

	print "Read OK"
	return N, den

########################################################################

def genTag(folder, step):
	filename = step.part_path
	outname = step.halo_path
	try :
		file = open(outname + ".tag", "rb")
		file.close()
	except  IOError:
		commande =  "./" + folder + "regroup -root %s -douter 80. -dsaddle 200. -dpeak 240. -f77 -o %s"%(outname,outname)
		print commande
		os.system(commande)
		
def readTag(filename):

	tagname = filename[:-10]  + "halo" + filename[-6:] + ".tag"

	file = open(tagname, "rb")
	print "Reading file", tagname

	dummy 	= np.fromfile(file, dtype=np.int32   ,count=1)	
	npart 	= np.fromfile(file, dtype=np.int32   ,count=1)[0]	
	ngrp 	= np.fromfile(file, dtype=np.int32   ,count=1)[0]
	dummy 	= np.fromfile(file, dtype=np.int32   ,count=1)
	dummy 	= np.fromfile(file, dtype=np.int32   ,count=1)
	tag 	= np.fromfile(file, dtype=np.int32   ,count=npart)
	dummy 	= np.fromfile(file, dtype=np.int32   ,count=1)

	file.close()
	
	print "Read OK"
	return npart, ngrp, tag
	
	
########################################################################

def genPos(folder, outname):
	try :
		file = open(outname + ".pos", "rb")
		file.close()
	except  IOError:
		commande =  "./" + folder + "poshalo -inp %s -pre zregroup -xmi 0.1 -xma 0.9"%(outname)
		print commande
		os.system(commande)

def readPos(filename):
#TODO
	return 0
