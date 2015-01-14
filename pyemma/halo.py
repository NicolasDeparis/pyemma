import os,sys
import numpy as np
import time
import matplotlib.pylab as plt

import IO
import part
import kdtree as kd
import hop
import physique


def readAll(filename):

	npartden, den       = hop.readDen(filename)
	nparttag, ngrp, tag = hop.readTag(filename)
	npart,a,par        = part.read(filename)
	nstar,a,star        = part.read(filename)

	return den, tag, par, star, nstar, a,ngrp

def getCenter(halo, den):

	mask = np.argmax(den)
	x  = halo.x[mask]
	y  = halo.y[mask]
	z  = halo.z[mask]
	
	return np.array((x,y,z), dtype = np.float32)
	
def getRvir(halo):		

	M = np.sum(halo.mass)
	return np.power(3.*M/(800*np.pi),1./3.)

def getHalo(par, tag, den, HALO):

	mask  = np.where(tag == HALO)
	nmask = len(mask[0])

	halo = part.Part(nmask,0)
	halo.mask(par,mask)

	return halo, getCenter(halo, den[mask]), getRvir(halo)
	
def findStar(center,Rvir, star, tree):
 
 	if tree : 
		mask3,dump = tree.locatenear(center,Rvir)			
	else : 
		mask3 = []
	
	nmask3 = len(mask3)
	
	starshalo = part.Part(nmask3,1)

	if nmask3 != 0:
		print nmask3
		starshalo.mask(star,mask3)

	return starshalo, np.array(mask3, dtype = np.int32)
	
def genHaloFiles(file):

	try :
		outname = file[:-10]  + "halo" + file[-6:] + ".halo"
		f = open(outname,'rb')	
		f.close()

	except  IOError:

		t0 = time.time()	
		den, tag, part, star, nstar, a, ngrp = readAll(file)	

		t1 = time.time()
		if nstar :	
			print "Tree generation"
			tree = kd.Tree(star.x,star.y,star.z)
			print "Tree generation OK"
		else :
			tree = 0

		t2 = time.time()
		print "Computation"

		center   = np.empty((ngrp), dtype=np.object)
		mask     = np.empty((ngrp), dtype=np.object)

		nstarhalo= np.zeros(ngrp, dtype=np.int32)
		Rvir     = np.zeros(ngrp, dtype=np.float32)
		Mh       = np.zeros(ngrp, dtype=np.float32)
		Ms       = np.zeros(ngrp, dtype=np.float32)


		for HALO in range(ngrp):

			halo,center[HALO],Rvir[HALO]   = getHalo(part,tag,den,HALO)
			starHalo, mask[HALO]           = findStar(center[HALO],Rvir[HALO],star, tree)

			nstarhalo[HALO] = len(mask[HALO])
			Mh[HALO] = np.sum(halo.mass)
			Ms[HALO] = np.sum(starHalo.mass)

		outname = file[:-10]  + "halo" + file[-6:] + ".halo"
		f = open(outname,'wb')		

		ngrp.tofile(f)		

		Rvir.tofile(f)
		Mh.tofile(f)
		Ms.tofile(f)
		nstarhalo.tofile(f)
		for i in range(ngrp):
			center[i].tofile(f)
			mask[i].tofile(f)			
		f.close()

		t3 = time.time()	

		print "Computation OK"
		print "Lecture", t1-t0
		print "tree",    t2-t1	
		print "Calcul",  t3-t2
		
		
def readHaloFile(file):

	outname = file[:-10]  + "halo" + file[-6:] + ".halo"
	try :
		f = open(outname,'rb')	
	except  IOError:
		print outname, "not found"

	ngrp = np.fromfile(f, dtype = np.int32, count = 1)		
	print ngrp

	center   = np.empty((ngrp), dtype=np.object)
	Rvir     = np.zeros( ngrp)
	Mh       = np.zeros( ngrp)
	Ms       = np.zeros( ngrp)	
	N       = np.zeros( ngrp)	
	mask     = np.empty((ngrp), dtype=np.object)


	Rvir   = np.fromfile(f, dtype = np.float32, count =ngrp)
	Mh     = np.fromfile(f, dtype = np.float32, count =ngrp)
	Ms     = np.fromfile(f, dtype = np.float32, count =ngrp)
	N      = np.fromfile(f, dtype = np.int32,   count =ngrp)

	for i in range(ngrp):
	
		center[i] = np.fromfile(f, dtype = np.float32, count = 3)
		mask[i]   = np.fromfile(f, dtype = np.int32, count =N[i])			

	f.close()


	return center,Rvir,Mh,Ms,N,mask

def Ms_f_Mh(Mh0,Ms0):

		param =  "../krusty/Quartz/data/wSN7/"
		Mh = physique.m2mo(Mh0,param)
		Ms = physique.m2mo(Ms0,param)

		print Mh
		print Ms

		b = 16
		n0, bins0 = np.histogram(Mh,bins=b)
		n1, bins1 = np.histogram(Ms,bins=b)
		

		x = np.zeros(b)
		y = np.zeros(b)
		for i in range(b):
			x[i] = n0[i] * bins0[i] #+ ( bins0[i+1] - bins0[i] )/2.
			y[i] = n1[i] * bins1[i] #+ ( bins1[i+1] - bins1[i] )/2.

		print x,y
		plt.clf()
		plt.loglog(Mh,Ms,'.')
		plt.title("Stars mass function of halo mass ")
		plt.xlabel(r'halo mass (M0)')
		plt.ylabel(r'stars mass (M0)')

		plt.legend()
		plt.show(block=False)



def plotHalo(file):
		
	fig = plt.figure()	
	
	ax = fig.add_subplot(1, 1, 1)

	den, tag, part, star, nstar, a, ngrp = readAll(file)	

	#tree = kd.Tree(star.x,star.y,star.z)

	center   = np.empty((ngrp), dtype=np.object)
	Rvir     = np.zeros(ngrp, dtype=np.float32)
	mask     = np.empty((ngrp), dtype=np.object)


	for HALO in range(ngrp):

		halo,center[HALO],Rvir[HALO]   = getHalo(part,tag,den,HALO)
		#starHalo, mask[HALO]           = findStar(center[HALO],Rvir[HALO],star, tree)

		#plt.plot(starHalo.x[HALO],starHalo.y[HALO], '.')
		
		ax.add_patch( plt.Circle((center[HALO][0],center[HALO][1]),Rvir[HALO],color='b',fill=False))
	#	plt.plot(center[HALO][0],center[HALO][1], 'ro')

	plt.show(block=False)
	
def getResolution(file):
	
	folder,filename = IO.splitPath(file)
	folder +="/"
	
	nstar,a,star = part.read(file)
	print "number of stars : %s"%str(nstar)
	L=10
	resmstarmin = physique.m2mo(np.min(star.mass),folder)
	resmstarmax = physique.m2mo(np.min(star.mass),folder)
	print "mass star = %s -> %s "%(str(resmstarmin),str(resmstarmin))

	file=file.replace("star","part")
	npart,a,par = part.read(file)
	print "number of DM part : %s"%str(npart)
	L=10
	resmpar = physique.m2mo(np.min(par.mass),folder)
	print "mass of DM star = %s"%str(resmpar)



