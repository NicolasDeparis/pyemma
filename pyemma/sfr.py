import numpy as np
import matplotlib.pylab as plt

import physique
import part
import IO
import halo
import kdtree as kd


def getFromPop(stars,folder,weight):
	b=32
	z=np.zeros(b)
	sfr=np.zeros(b)
	
	
	if len(stars.mass):
		param = IO.ParamsInfo(folder = folder).get()
		L = float(param["unit_l"])/3.085677e16 /1e6	

		n0,bin_edges=np.histogram(stars.age,bins=b)
	
		z=physique.a2z(physique.t2a( bin_edges))[:-1]

		n=len(n0)
		for i in range(n):
			sfr[i] = physique.m2mo( stars.mass[i] * (n0[i])/weight ,folder )  / float( bin_edges[i+1] - bin_edges[i]) / pow(L,3)
	
	return z,sfr
	
def fromPop(stars,folder,lab,weight):
		z,sfr = getFromPop(stars,folder,weight)	
		if len(stars.mass):
			plt.semilogy(z,sfr, label = lab)
			plt.xlabel('z')
			plt.ylabel(r'$M_{\odot}.yrs^{-1}.(Mpc.h^{-1})^{-3}$')

def fromSnap(file, lab):
	folder,filename = IO.splitPath(file)
	Ntot,a,stars=part.read(file)
	fromPop(stars,folder,lab,1)
	plt.legend()
	plt.show(block=False)

def observation():
	"data from Bouwens et al. 2008"
	
	Z   = [0.2,2,3.8,4.9,5.9,6.8,7.9,10.4]

	SFR1 = [-2.2,-1.4,-1.44,-1.53,-1.84,-1.95,-2.2,-3.33]
	SFR1 = [ pow(10,x) for x in SFR1]
	plt.plot(Z,SFR1, 'k--')
#, label = "obs dust uncorrected"

	
	SFR2 = [-1.7,-1,-1.06,-1.19,-1.59,-1.72,-2.05,-3.18]
	SFR2 = [ pow(10,x) for x in SFR2]
	plt.plot(Z,SFR2, 'k--')
#, label = "obs dust corrected"

	plt.show(block=False)


def fromSnap_2(file, lab):

	folder,filename = IO.splitPath(file)

	if "star." in filename:
		star = 1
	else:
		star = 0
	s = 10 + star
	
	param = IO.Params(folder = folder).get()
	L = float(param["unit_l"])/3.085677e16 /1e6	

	nproc = IO.getNproc(file)
	
	for proc in range(nproc):
		Ntot,a,data =part.read1proc(file + ".p"+ str(proc).zfill(5))
		stars = part.Part(Ntot, star)
		i=0
		for j in range(0,data.shape[0],s):
			stars.define(data[j:j+s] ,i )
			i+=1
			
		if len(stars.mass):
			b=32
			n0,bin_edges=np.histogram(stars.age,bins=b)
		
			z=physique.a2z(physique.t2a( bin_edges))

			n=len(n0)
			sfr=np.zeros(n)
			for i in range(n):
				sfr[i] = physique.m2mo( stars.mass[0] * (n0[i]) ,folder )  / float( bin_edges[i+1] - bin_edges[i]) / pow(L,3)
				
			plt.semilogy(z[:-1],sfr)
			
			plt.xlabel('z')
			plt.ylabel(r'$M_{\odot}.yrs^{-1}.(Mpc.h^{-1})^{-3}$')
			plt.ylim(1e-6,1)
			
			plt.legend()
	plt.show(block=False)
	
	
	
def haloMass(file):
	folder,filename = IO.splitPath(file)

	den, tag, part, star, nstar, a, ngrp = halo.readAll(file)	

	if nstar :	
		print "Tree generation"
		tree = kd.Tree(star.x,star.y,star.z)
		print "Tree generation OK"
	else :
		tree = 0

	print "Computation"

	center   = np.empty((ngrp), dtype=np.object)
	mask     = np.empty((ngrp), dtype=np.object)
	starHalo = np.empty((ngrp), dtype=np.object)
	
	nstarhalo= np.zeros(ngrp, dtype=np.int32)
	Rvir     = np.zeros(ngrp, dtype=np.float32)
	Mh       = np.zeros(ngrp, dtype=np.float32)



	for HALO in range(ngrp):
		h,center[HALO],Rvir[HALO]= halo.getHalo(part,tag,den,HALO)
		starHalo[HALO], mask[HALO]= halo.findStar(center[HALO],Rvir[HALO],star, tree)

		Mh[HALO] = np.sum(h.mass)


	return Mh, starHalo
	
def haloMassPlot(Mh, starHalo, folder):

	Nhalo=len(Mh)
	
	Nbins=8
	n,bin_edges=np.histogram(Mh,bins=Nbins)

	sfr = np.empty(Nbins, dtype=np.object)
	
	for i in range(Nbins):
		sfr[i] = np.zeros(32)	
	
	for b in range(Nbins):
		star= part.Part(0,1)
		weight=0
		for HALO in range(Nhalo):
			if Mh[HALO] > bin_edges[b] and Mh[HALO] < bin_edges[b+1]:
				star.append(starHalo[HALO])
				weight +=1
	
		lab = "%s<Mhalo<%s"%("{:.2e}".format(physique.m2mo(bin_edges[b], folder)), "{:.2e}".format(physique.m2mo(bin_edges[b+1],folder)))
		fromPop(star,folder,lab,weight)
		
	plt.legend()
	plt.show(block=False)
