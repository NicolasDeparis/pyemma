import numpy as np
import os

import IO

def genHopFiles(filename):

	if not "part." in filename:
		print "HOP need a particle file"
		sys.exit(1)

	nproc = IO.getNproc(filename)
	folder  = "utils/hop/"

	outname = filename[:-10]  + "halo" + filename[-6:]

	try :
		file = open(outname + ".den", "rb")
		file.close()
	except  IOError:
		commande="./" + folder + "hop -in " + filename + " -o " + outname +" -p 1 -nf " + str(nproc)
		print commande
		os.system(commande)


	try :
		file = open(outname + ".tag", "rb")
		file.close()
	except  IOError:
		#-douter 80. -dsaddle 200. -dpeak 240.
		#-douter 20. -dsaddle 35. -dpeak 50.
		commande =  "./" + folder + "regroup -root %s -douter 80. -dsaddle 200. -dpeak 240. -f77 -o %s"%(outname,outname)
		print commande
		os.system(commande)

def readDen(filename):

	denname = filename[:-10]  + "halo" + filename[-6:] + ".den"

	file = open(denname, "rb")
	print "Reading file", denname

	N = np.fromfile(file, dtype=np.int32   ,count=1)[0]
	den = np.fromfile(file, dtype=np.float32 ,count=N)
	
	file.close()

	print "Read OK"
	return N, den

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
