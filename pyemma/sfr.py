import numpy as np
import matplotlib.pylab as plt

import physique
import part
import IO
import halo
import kdtree as kd
import param

def getFromPop(stars,folder):
	
	Nbins=8

	z=np.zeros(Nbins)
	sfr=np.zeros(Nbins)		

	p = param.ParamsInfo(folder = folder).get()
	L = float(p["unit_l"])/3.085677e16 /1e6	
	dv = pow(L,3)
	
	age_min = np.min(stars.age)
	age_max = np.max(stars.age)
	
	bins_edges = np.linspace(age_min,age_max,Nbins+1)
	weights=physique.m2mo(stars.mass,folder)
	
	sfr,bin_edges=np.histogram(stars.age,bins=bins_edges, weights=weights)
	sfr/=np.diff(bin_edges)
	
	z_all=physique.a2z(physique.t2a(bin_edges))
	z=(z_all[1:]+z_all[:-1])/2		

	return z,sfr
	
def fromSnap_cpu(file, lab):

	folder,filename = IO.splitPath(file)
	nproc = IO.getNproc(file)

	param = param.ParamsInfo(folder = folder).get()
	L = float(param["unit_l"])/3.085677e16 /1e6	
	dv = pow(L,3)
	
	for proc in range(nproc):
		print "cpu=%d"%proc
		N,a,data =part.read1proc(file + ".p"+ str(proc).zfill(5))
		print N
		stars = part.Part(N, 1)
		
		ii=0
		for j in range(0,data.shape[0],11):
			stars.define(data[j:j+11] ,ii )
		
		if N:
			z,sfr=getFromPop(stars,folder,1)
	
			plt.semilogy(z,sfr)		
			plt.plot(z,sfr)
			plt.xlabel('z')
			plt.ylabel(r'$M_{\odot}.yrs^{-1}.(Mpc.h^{-1})^{-3}$')		
			plt.legend()
					
	plt.show(block=False)
	
def fromPop(stars,folder,**kwargs):
	z,sfr = getFromPop(stars,folder)	
	if len(stars.mass):
		plt.semilogy(z,sfr,**kwargs)
		#plt.plot(z,sfr,**kwargs)
		plt.xlabel('z')
		plt.ylabel(r'$M_{\odot}.yrs^{-1}.(Mpc.h^{-1})^{-3}$')
		

def fromSnap(file, **kwargs):
#	plt.figure()
	folder,filename = IO.splitPath(file)
	Ntot,a,stars=part.read(file)
	fromPop(stars,folder,**kwargs)
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


def haloMass(file):
	folder,filename = IO.splitPath(file)
	den, tag, part, star, nstar, a, ngrp = halo.readAll(file)

	tree = halo.genTree(star)

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
#histogtram sur les masses de halo
	
def haloMassPlot(Mh, starHalo,folder,**kwargs):

	cut1=physique.mo2m(5e9,folder)
	cut2=physique.mo2m(1e11,folder)
	bin_edges=[0,cut1,cut2]
	Nbins=len(bin_edges)-1

	sfr = np.empty(Nbins, dtype=np.object)

	for i in range(Nbins):
		sfr[i] = np.zeros(32)
	
	for b in range(Nbins):
		star= part.Part(0,1)
		weight=0
		for HALO in range(len(Mh)):
			if Mh[HALO] > bin_edges[b] and Mh[HALO] < bin_edges[b+1]:
				star.append(starHalo[HALO])
				weight +=1
	
		mmin=physique.m2mo(bin_edges[b], folder)		
		mmax=physique.m2mo(bin_edges[b+1], folder)
		labmin="{:.0e}".format(mmin)
		labmax="{:.0e}".format(mmax)
		kwargs["label"]="%s<Mh<%s N=%s"%(labmin,labmax,weight)	
		fromPop(star,folder,**kwargs)
	plt.xlim(6,13)
	plt.legend()
	plt.show(block=False)

