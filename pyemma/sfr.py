# coding: utf8

import numpy as np
import matplotlib.pylab as plt

import physique
import part
import IO
import halo
import kdtree as kd
import param

def getFromPop(stars,Param,Nbins=10):
	
	z=np.zeros(Nbins)
	sfr=np.zeros(Nbins)		

	p = Param.info.get()
	L = float(p["unit_l"])/3.085677e16/1e6
	dv = pow(L,3)
	
	age_min = np.min(stars.age)
	age_max = np.max(stars.age)
	dt = age_max - age_min
	bins_edges = np.linspace(age_min,age_max,Nbins+1)
	weights=physique.m2mo(stars.mass,Param)
	
	mass = np.unique(weights)
	if len(mass)>2:
		print "WARNING : there an ambiquity about the masses of stars"
		print "mass = ", mass
	
	mass = np.max(mass)
	#mass *= 10
	weights =0*weights+mass	
	
	sfr,bin_edges=np.histogram(stars.age, bins=bins_edges, weights=weights)
	sfr/=np.diff(bin_edges)*dv

	z_all=physique.a2z(physique.t2a(bin_edges))
	z=(z_all[1:]+z_all[:-1])/2		
	
	avg_sfr=len(stars.age)*mass/dt/dv
	print "average SFR", avg_sfr
	return z,sfr
	

	
def fromSnap_cpu(step, lab):

#	folder,filename = IO.splitPath(file)

	nproc = IO.getNproc(step.star_path)

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
			plt.ylabel(r'$M_{\odot}.yrs^{-1}.Mpc^{-3}$')		
			plt.legend()
					
	plt.show(block=False)
	
def fromPop(stars,Param,**kwargs):
	z,sfr = getFromPop(stars,Param)	
	if len(stars.mass):
		#plt.semilogy(z,sfr,**kwargs)
		plt.plot(z,np.log10(sfr),**kwargs)
		#plt.plot(z,sfr,**kwargs)
		plt.xlabel('z')
		plt.ylabel(r'$ LOG_{10} (SFR) [ M_{\odot}.yrs^{-1}.Mpc^{-3} ]$')		
			

def fromSnap(step, **kwargs):
#	plt.figure()

	Ntot,a,stars=part.read(step.star_path)
	P=param.Param(step.step_folder,step.number)
	
	fromPop(stars,P,**kwargs)
	plt.legend()
	plt.show(block=False)

def from_info(folder):
	dat = param.ParamsAvg(0,folder).get()

	print dat

def observation1():
	"""
	"data from Bouwens et al. 2008"
	#Z   = [0.2,2,3.8,4.9,5.9,6.8,7.9,10.4]
	#SFR1 = [-2.2,-1.4,-1.44,-1.53,-1.84,-1.95,-2.2,-3.33]
	Z   = [0.2,2,3.8,4.9,5.9,6.8,7.9]
	SFR1 = [-2.2,-1.4,-1.44,-1.53,-1.84,-1.95,-2.2]
	#SFR1 = [ pow(10,x) for x in SFR1]
	plt.plot(Z,SFR1, 'ro', label = "Bouwens et al. 2015")
#, label = "obs dust uncorrected"
	#SFR2 = [-1.7,-1,-1.06,-1.19,-1.59,-1.72,-2.05,-3.18]
	SFR2 = [-1.7,-1,-1.06,-1.19,-1.59,-1.72,-2.05]
	#SFR2 = [ pow(10,x) for x in SFR2]
	plt.plot(Z,SFR2, 'ro')
#, label = "obs dust corrected"
	plt.show(block=False)
	plt.legend()
	"""
	
	"""
	data from Bouwens et al. 2015 http://arxiv.org/pdf/1403.4295v4.pdf

	format: 
		z, log10(sfr)_dust_uncorrected, error_inf, error_sup, log10(sfr)_dust_corrected, error_inf, error_sup
	"""	
	z=np.array([
	[3.8,	-1.38, 	0.06, 0.06,	-1.00, 	0.06,	0.06],
	[4.9,  	-1.60,	0.06, 0.06,	-1.26, 	0.06,	0.06],
	[5.9, 	-1.80,	0.06, 0.06,	-1.55, 	0.06,	0.06],
	[6.8,  	-1.92,	0.06, 0.06,	-1.69, 	0.06,	0.06],
	[7.9,  	-2.23,	0.07, 0.06,	-2.08, 	0.07,	0.07],
	[10.4, 	-3.28,	0.36, 0.45, -3.13,	0.36,	0.45]])
	

	x=z[:,0]
	
	y1=np.power(10., z[:,1])
	low_err1 = np.abs(y1 - np.power(10., z[:,1] - z[:,3]))
	sup_err1 = np.abs(y1 - np.power(10., z[:,1] + z[:,2]))
	yerror1 = [low_err1, sup_err1]
	
	


	y2=np.power(10., z[:,4])
	low_err2 = np.abs(y2 - np.power(10., z[:,1] - z[:,6]))
	sup_err2 = np.abs(y2 - np.power(10., z[:,1] + z[:,5]))
	yerror2 = [low_err2, sup_err2]
	
	#plt.errorbar(x, y, yerr=yerror, ls='none', fmt='b*',label = "Bouwens et al. 2015")	
	#plt.errorbar(x,y, yerr=yerror, ls='none', fmt='b*')	

	x = np.append(x,x)
	y = np.append(y1,y2)
	
	low_err = np.append(low_err1,low_err2)
	sup_err = np.append(sup_err1,sup_err2)
	yerror = [low_err, sup_err]
	
	
	print x,y, yerror

	return x, y, None, yerror
	
	
	

def observation2():
	"""
	Observational constraints for SFR.
	data from  Madau, P., & Dickinson, M. 2014, ARA&A, 52, 415 (MD 14)

	format: 
		zmin, zmax, log10(sfr), error_inf, error_sup
	"""
	z=np.array([
	[0.01,	0.1, 	-1.82,	+0.09,	-0.02],
	[0.2,	0.4,  	-1.50,	+0.05,	-0.05],
	[0.4,	0.6,  	-1.39,	+0.15,	-0.08],
	[0.6,	0.8,  	-1.20,	+0.31,	-0.13],
	[0.8,	1.2,  	-1.25,	+0.31,	-0.13],
	[0.05,	0.05,	-1.77,	+0.08,	-0.09],
	[0.05,	0.2,  	-1.75,	+0.18,	-0.18],
	[0.2,	0.4,  	-1.55,	+0.12,	-0.12],
	[0.4,	0.6,  	-1.44,	+0.10,	-0.10],
	[0.6,	0.8,  	-1.24,	+0.10,	-0.10],
	[0.8,	1.0,  	-0.99,	+0.09,	-0.08],
	[1.0,	1.2,  	-0.94,	+0.09,	-0.09],
	[1.2,	1.7,  	-0.95,	+0.15,	-0.08],
	[1.7,	2.5,  	-0.75,	+0.49,	-0.09],
	[2.5,	3.5,  	-1.04,	+0.26,	-0.15],
	[3.5,	4.5,  	-1.69,	+0.22,	-0.32],
	[0.92,	1.33, 	-1.02,	+0.08,	-0.08],
	[1.62,	1.88, 	-0.75,	+0.12,	-0.12],
	[2.08,	2.37, 	-0.87,	+0.09,	-0.09],
	[1.9,	2.7, 	-0.75,	+0.09,	-0.11],
	[2.7,	3.4, 	-0.97,	+0.11,	-0.15],
	[3.8,	3.8, 	-1.29,	+0.05,	-0.05],
	[4.9,	4.9, 	-1.42,	+0.06,	-0.06],
	[5.9,	5.9, 	-1.65,	+0.08,	-0.08],
	[7.0,	7.0, 	-1.79,	+0.10,	-0.10],
	[7.9,	7.9, 	-2.09,	+0.11,	-0.11],
	[7.0,	7.0, 	-2.00,	+0.10,	-0.11],
	[8.0,	8.0, 	-2.21,	+0.14,	-0.14],
	[0.03, 	0.03, 	-1.72,	+0.02,	-0.03],
	[0.03,	0.03, 	-1.95,	+0.20,	-0.20],
	[0.40,	0.70,  	-1.34,	+0.22,	-0.11],
	[0.70,	1.00,  	-0.96,	+0.15,	-0.19],
	[1.00,	1.30,  	-0.89,	+0.27,	-0.21],
	[1.30,	1.80,  	-0.91,	+0.17,	-0.21],
	[1.80,	2.30,  	-0.89,	+0.21,	-0.25],
	[0.40,	0.70,  	-1.22,	+0.08,	-0.11],
	[0.70,	1.00,  	-1.10,	+0.10,	-0.13],
	[1.00,	1.30,  	-0.96,	+0.13,	-0.20],
	[1.30,	1.80,  	-0.94,	+0.13,	-0.18],
	[1.80,	2.30,  	-0.80,	+0.18,	-0.15],
	[0.00,	0.30,  	-1.64,	+0.09,	-0.11],
	[0.30,	0.45,  	-1.42,	+0.03,	-0.04],
	[0.45,	0.60,  	-1.32,	+0.05,	-0.05],
	[0.60,	0.80,  	-1.14,	+0.06,	-0.06],
	[0.80,	1.00,  	-0.94,	+0.05,	-0.06],
	[1.00,	1.20,  	-0.81,	+0.04,	-0.05],
	[1.20,	1.70,  	-0.84,	+0.04,	-0.04],
	[1.70,	2.00,  	-0.86,	+0.02,	-0.03],
	[2.00,	2.50,  	-0.91,	+0.09,	-0.12],
	[2.50,	3.00,  	-0.86,	+0.15,	-0.23],
	[3.00,	4.20,  	-1.36,	+0.23,	-0.50],
	])
	
	
	x= (z[:,1]+z[:,0])/2.
	xerror = (z[:,1]-z[:,0])/2.
	
	y= np.power(10., z[:,2])
		
	low_err = np.abs(y - np.power(10., z[:,2] + z[:,4]))
	sup_err = np.abs(y - np.power(10., z[:,2] + z[:,3]))	
	yerror = [low_err, sup_err]
	
	return x, y, xerror, yerror
	
	#plt.errorbar(x, y, xerr=xerror, yerr=yerror, ls='none', fmt='b8', label = "Madau et al. 2014")	
	
	
def observation():
	x1, y1, xerror1, yerror1 = observation1()
	x2, y2, xerror2, yerror2 = observation2()

	separe=0
	if separe:
		plt.errorbar(x1, y1, xerr=xerror1, yerr=yerror1, ls='none', fmt='b8', label = "Madau et al. 2014")	
		plt.errorbar(x2, y2, xerr=xerror2, yerr=yerror2, ls='none', fmt='b8', label = "Bouwens et al. 2015")	
	else:
		plt.errorbar(x1, y1, xerr=xerror1, yerr=yerror1, ls='none', fmt='k.', label = "Observation")	
		plt.errorbar(x2, y2, xerr=xerror2, yerr=yerror2, ls='none', fmt='k.')

def xion_constrain():
	"""
	Observational constraints for xion.
	data from Bouwens 2015
	"""
	
	
	#er=[[0.5],[1]]
	#plt.errorbar(6., 0.01, xerr=er, ls='none', fmt='ko')
	#plt.errorbar(6., 0.5,  xerr=er, ls='none', fmt='ko')
	
	z=np.array([

	[5.03, 0.9999451, 0.0000142, -0.0000165],
	[5.25, 0.9999330, 0.0000207, -0.0000244],
	[5.45, 0.9999333, 0.0000247, -0.0000301],
	[5.65, 0.9999140, 0.0000365, -0.0000460],
	[5.85, 0.9998800, 0.0000408, -0.0000490],
	[6.10, 0.99957  , 0.00030  , -0.00030  ],
	[7.0 , 0.66     , 0.12     , -0.09 	   ]
	])
	
	x=z[:,0]
	y=z[:,1]

	yerror = [-z[:,3], z[:,2]]
	
	
	plt.errorbar(x, 1.-y, xerr=None, yerr=yerror, ls='none', fmt='ko', label="Observations")		
	
	
	"""
	5.9 >0.89 Dark Gaps in Quasar Spectra McGreer et al. (2015)
	5.6 >0.91 Dark Gaps in Quasar Spectra McGreer et al. (2015)
	6.24-6.42 <0.9 (2) Ly Damping Wing of Quasars Schroeder et al. (2013)
	Higher-Redshift Constraints
	
	
	Prevalence of Lyα Emission in Galaxies S14
	4. 8.0 QHII (z = 8) < 0.35 Prevalence of Lyα Emission in Galaxies S14
	Continuity with Ionizing Emissivity Estimates at z = 4.75
	5. log10 N˙
	ion(z = 4.75) = 1050.99±0.45 s
	−1 Mpc−3 BB13
	Other Constraints on the Ionization History of the Universe Not Explicitly Usedb
	
	6.3 ≥0.5 Lyα Damping Wing of a GRB Totani et al. (2006)	McQuinn et al. (2008)
	6.6 ≥0.6 Lyα Emitters Ouchi et al. (2010)
	6.6 ≥0.5 Galaxy Clustering McQuinn et al. (2007),	Ouchi et al. (2010)
	
	7.0 0.32-0.64 Lyα-Emitter LFs Ota et al. (2008)
	7.0 ∼0.5 Prevalence of Lyα Emission in Galaxies Caruana et al. (2014)
	7.0 0.1-0.4 Prevalence of Lyα Emission in Galaxies Ono et al. (2012)
	7.0 <0.49 Prevalence of Lyα Emission in Galaxies P14
	7.0 <0.5 Prevalence of Lyα Emission in Galaxies R13c
	7.0 <0.5 Clustering of Lyα Emitting Galaxies Sobacchi & Mesinger (2015)
	7.1 ≤0.9 Near-Zone Quasar Mortlock et al. (2011),
	Bolton et al. (2011)
	8.0 <0.70 Prevalence of Lyα Emission in Galaxies Tilvi et al. (2014)
	a This table is a compilation of the constraints presented in the original papers under References, but with valuable guidance
	"""
	
	
"""
#if z < zp 1 − QHII(z) ∝ (1 + z)
#if z ≥ zp QHII(z) ∝ exp(−λ(1 + z)). (2)

zp=6.1
l=0.73
qzp = 0.99986



N=100
Q=np.zeros(N)
i=0

z = np.linspace(4,14,N)
for curz in z:
    if curz<zp:
        Q[i]=1+ (1+curz)**3
    else:
        Q[i]=np.exp(-l*(1+curz))*30
        
    i+=1
plt.plot(z,Q)
"""   
	
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
#histogtramme sur les masses de halo
	
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

def roberson_2015_fit():
    """ from : Cosmic Reionization and Early Star-forming Galaxies: A Joint Analysis of New Constraints from Planck and the Hubble Space Telescope"""
    z=np.linspace(0,20,1000)
    
    """Plank case"""
    a=0.01376
    b=3.26
    c=2.59
    d=5.68
    y_fit= a*np.power(1+z, b)/(1+np.power((1+z)/c,d))
    plt.plot(z,y_fit)
    
    """WMAP case"""
    a=0.01306
    b=3.66
    c=2.28
    d=5.29
    y_fit= a*np.power(1+z, b)/(1+np.power((1+z)/c,d))
    #plt.plot(z,y_fit)
       
    
def roberson_2010_fit():
    """from : Early star-forming galaxies and the reionization of the Universe"""
    z=np.linspace(0,20,100)
    
    a = 0.009
    b = 0.27
    h = 2.5
    
    """metal poor case"""
    c = 3.7
    d = 7.4
    g = 1e-3
    y_fit  = (a + b*(z/c)**h ) /(1. + (z/c)**d ) + g    
    plt.plot(z,y_fit)

    """metal rich case"""
    c = 3.4
    d = 8.3
    g = 1e-4
    y_fit  = (a + b*(z/c)**h ) /(1. + (z/c)**d ) + g
    plt.plot(z,y_fit)
