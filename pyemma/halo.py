import os,sys
import numpy as np
import time
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle

import IO
import part
import kdtree as kd
import hop
import physique
import param
import part
import plot as p

def getStarburstNphot():
	data = np.loadtxt('python/Starburst99/912_inst_e.dat.txt',unpack=True)
	x=data[0]
	y=np.power(10,data[1])
	return x,y

class Halo:
	def __init__(self,ngrp):
		self.t=0
		self.file=""
		
		self.N=ngrp
		self.x=np.zeros(ngrp, dtype=np.float32)
		self.y=np.zeros(ngrp, dtype=np.float32)
		self.z=np.zeros(ngrp, dtype=np.float32)
		self.R=np.zeros(ngrp, dtype=np.float32)
		self.part=np.empty(ngrp, dtype=np.object)
		self.star=np.empty(ngrp, dtype=np.object)

	def getN(self):
		return self.N
	def toFile(self,name):
		with open(name, 'wb') as output:		
			pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)			

	def plotNumber(self,ax,n,axe1="x",axe2="y",**kwargs):	
		axes = {}
		axes["x"] = self.x[n]
		axes["y"] = self.y[n]
		axes["z"] = self.z[n]
		x=axes[axe1]+ self.R[n]
		y=axes[axe2]
		bbox_args = dict(boxstyle="round", fc="1")
		ax.annotate("%d"%n, xy=(x,y),bbox=bbox_args)
		
	def plot(self,ax,n,axe1="x",axe2="y",**kwargs):
		axes = {}
		axes["x"] = self.x[n]
		axes["y"] = self.y[n]
		axes["z"] = self.z[n]
		x=axes[axe1]
		y=axes[axe2]
		ax.add_patch(plt.Circle((x,y),self.R[n],fill=False, **kwargs))
	
	def plotPart(self,n,axe1="x",axe2="y",**kwargs):
		axes = {}
		axes["x"] = self.part[n].x
		axes["y"] = self.part[n].y
		axes["z"] = self.part[n].z
		x=axes[axe1]
		y=axes[axe2]
		plt.plot(x,y,'.',**kwargs)
		
	def plotStar(self,n,axe1="x",axe2="y",**kwargs):
		axes = {}
		axes["x"] = self.star[n].x
		axes["y"] = self.star[n].y
		axes["z"] = self.star[n].z
		x=axes[axe1]
		y=axes[axe2]
		plt.plot(x,y,'.',**kwargs)
		
	def plotPart3D(self,ax,n,**kwargs):
		ax.scatter(self.part[n].x,self.part[n].y,self.part[n].z,'.',**kwargs)
	def plotStar3D(self,ax,n,**kwargs):
		ax.scatter(self.star[n].x,self.star[n].y,self.star[n].z,'.',**kwargs)

	def getNphot(self,n,t0):
		N_dot=0
		t,Nphot = getStarburstNphot()
		
		for age in self.star[n].age:				
			T=t0-age
			photStar=Nphot[(np.abs(t-T)).argmin()]/1e50
			N_dot+=photStar
		return N_dot*1e50


	
def readHalo(name):	
	with open(name, 'rb') as input:
		return pickle.load(input)			
		
def getCenter(halo,den):
	"""return the center of a halo"""
	pos=np.argmax(den)
	x  = halo.x[pos]
	y  = halo.y[pos]
	z  = halo.z[pos]
	return np.array((x,y,z), dtype = np.float32)
	
def getR200(halo):
	"""return the R200 of a halo"""
	M = np.sum(halo.mass)
	return np.power(3.*M/(200.*4.*np.pi),1./3.)*5

def getHalo(par,mask):
	"""return all the particle of a halo"""
	halo = part.Part(len(mask[0]),0)
	halo.mask(par,mask)
	return halo

def getStar(center,R200,star,tree):
	"""find star belong to a halo"""
	mask,dump = tree.locatenear(center,R200)
	stars_halo = part.Part(len(mask),1)
	stars_halo.mask(star,mask)
	return stars_halo

def genTree(file,part,type):
	"""compute kd tree if it doensn't already exist, else read it"""
	name=file.replace("part.","halo.")+ ".%stree"%type
	try :
		open(name)
		tree=kd.Tree.FromFile(name)
	except  IOError:
		tree=kd.Tree(part.x,part.y,part.z)
		tree.WriteToFile(name)	
	return tree

def genHaloFiles(file):
	"""reduce data and save file"""
	name=file.replace("part.","halo.")+ ".halo"
	try :
		with open(name): pass
	except  IOError:

		npartden, den = hop.readDen(file)
		nparttag, ngrp, tag = hop.readTag(file)
		npart,a,part_all = part.read(file)
		nstar,a,star_all = part.read(file.replace("part","star"))

		tree=genTree(file,star_all,"star")
		halo=Halo(ngrp)
	
		for iHalo in range(ngrp):
			mask  = np.where(tag==iHalo)
	
			part_halo=getHalo(part_all,mask)
			center = getCenter(part_halo,den[mask])
			R200 = getR200(part_halo)
			star_halo = getStar(center,R200,star_all,tree)

			halo.x[iHalo]=center[0]
			halo.y[iHalo]=center[1]
			halo.z[iHalo]=center[2]
			halo.R[iHalo]=R200
			halo.part[iHalo]=part_halo
			halo.star[iHalo]=star_halo
		halo.toFile(name)

####################################################################################################
####################################################################################################

def plotStar(star,axe1="x",axe2="y",**kwargs):
	axes = {}
	axes["x"] = star.x
	axes["y"] = star.y
	axes["z"] = star.z
	x=axes[axe1]
	y=axes[axe2]
	plt.plot(x,y,'.',**kwargs)

def plotAll(file,axe1="x",axe2="y",**kwargs):
	name=file.replace("part.","halo.")+ ".halo"
	halo = readHalo(name)
	nstar,a,star = part.read(file.replace("part","star"))
	
	mtot_star= np.sum(star.mass)

	fig=plt.figure()
	ax = fig.add_subplot(1,1,1)
	
	for iHalo in range(halo.getN()):
		halo.plot(ax,iHalo,axe1,axe2)
		halo.plotPart(iHalo,axe1,axe2,c='k')
	
	plotStar(star,axe1,axe2,c='r')
	
	mstar_halo=0
	for iHalo in range(halo.getN()):
		halo.plotStar(iHalo,axe1,axe2,c='g')
		halo.plotNumber(ax,iHalo,axe1,axe2,)
		mstar_halo+=np.sum(halo.star[iHalo].mass)	
		
	print nstar_halo/mtot_star
	plt.title("%s%s"%(axe1,axe2))
	plt.xlim(0,1)
	plt.ylim(0,1)


def plotOne(file,	halo_number,axe1="x",axe2="y", **kwargs):
	name=file.replace("part.","halo.")+ ".halo"
	halo = readHalo(name)

	nstar,a,star_all = part.read(file.replace("part","star"))
	tree=genTree(file,star_all,"star")
	
	center=(halo.x[halo_number],halo.y[halo_number],halo.z[halo_number])
	R200=halo.R[halo_number]*3
	
	star_halo = getStar(center,R200,star_all,tree)

	fig=plt.figure()
	ax = fig.add_subplot(1,1,1)

	halo.plot(ax,halo_number,axe1,axe2)
	halo.plotNumber(ax,halo_number,axe1,axe2,)
	
	halo.plotPart(halo_number,axe1,axe2,c='k')
	plotStar(star_halo,axe1,axe2,c='r')
	halo.plotStar(halo_number,axe1,axe2,c='g')
		

def compare(file1,file2,halo_number, **kwargs):
	
	name=file1.replace("part.","halo.")+ ".halo"
	halo1 = readHalo(name)
	
	name=file2.replace("part.","halo.")+ ".halo"
	halo2 = readHalo(name)
	
	tree=genTree(file1,halo1,"halo")

	x=halo2.x[halo_number]
	y=halo2.y[halo_number]
	z=halo2.z[halo_number]
	mask=tree.nearest((x,y,z))
	
	axe1="x"
	axe2="y"
	
	plotOne(file1,mask,axe1,axe2, **kwargs)
	plotOne(file2,halo_number,axe1,axe2, **kwargs)


def test3D(file, **kwargs):
	name=file.replace("part.","halo.")+ ".halo"
	halo = readHalo(name)
	
	halo_number=7

	fig=plt.figure()
	ax = fig.add_subplot(111,projection='3d')

	nstar,a,star = part.read(file.replace("part","star"))
	
	dx=halo.R[halo_number]*5
	x=halo.x[halo_number]
	y=halo.y[halo_number]
	z=halo.z[halo_number]
	
	condx=(star.x>x-dx)*(star.x<x+dx)
	condy=(star.y>y-dx)*(star.y<y+dx)
	condz=(star.z>z-dx)*(star.z<z+dx)
	
	cond=np.where(condx*condy*condz)
	
	ax.scatter(star.x[cond],star.y[cond],star.z[cond],'.',c='g')

	halo.plotPart3D(ax,halo_number,c='k')
	halo.plotStar3D(ax,halo_number,c='r')
	
	
		
########################################################################
## Luminosity fonction
########################################################################


def luminosity(file):
	folder,filename = IO.splitPath(file)
	p = param.ParamsInfo(folder).get()
	cst=physique.Constantes()
	u_mass=p["unit_mass"]/cst.M0
	
	name = np.where( [".p00000" in n for n in os.listdir(folder)])[0]
	a = part.geta(folder+os.listdir(folder)[name])
	t= physique.a2t(a)
	
	halo = readHalo(file.replace("part.","halo.")+ ".halo")
	n=halo.getN()
	
	mh=np.zeros(n, dtype=np.float64)
	nphot=np.zeros(n, dtype=np.float64)
	
	for iHalo in range(n):
		par=halo.part[iHalo]
		mh[iHalo]=np.sum(par.mass)*u_mass
		nphot[iHalo] = halo.getNphot(iHalo,t)
		print mh[iHalo], nphot[iHalo]
	
	
	Nbins=len(mh)
	bins_edges = np.linspace(np.min(mh),np.max(mh),Nbins+1)
	
	nphot,bin_edges=np.histogram(mh,bins=bins_edges,weights=nphot)
	mh=(bin_edges[1:]+bin_edges[:-1])/2

	
	print mh,nphot
	plt.loglog(mh,nphot,'o')
	
	
####################################################################################################
####################################################################################################



def Ms_f_Mh(Mh0,Ms0,folder,**kwargs):

		Mh = physique.m2mo(Mh0,folder)
		Ms = physique.m2mo(Ms0,folder)

		"""
		print Mh
		print Ms
		"""
		
		"""
		b = 2
		n0, bins0 = np.histogram(Mh,bins=b)
		n1, bins1 = np.histogram(Ms,bins=b)
		

		x = np.zeros(b)
		y = np.zeros(b)
		for i in range(b):
			x[i] = n0[i] * bins0[i] #+ ( bins0[i+1] - bins0[i] )/2.
			y[i] = n1[i] * bins1[i] #+ ( bins1[i+1] - bins1[i] )/2.
		"""
						
	#	plt.clf()
		plt.loglog(Mh,Ms,'.', **kwargs)
		plt.title("Stars mass function of halo mass")
		plt.xlabel(r'halo mass ($M_{\odot}$ )')
		plt.ylabel(r'stars mass ($M_{\odot}$ )')

		plt.legend()
		plt.show(block=False)



def plot(file, fig, folder,**kwargs):
	
	ax = fig.add_subplot(1, 1, 1)
	den, tag, part, star, nstar, a, ngrp = readAll(file)	

	center   = np.empty((ngrp), dtype=np.object)
	R200     = np.zeros(ngrp, dtype=np.float32)

	for HALO in range(ngrp):
		halo,center[HALO],R200[HALO]   = getHalo(part,tag,den,HALO)	
		if physique.m2mo(np.sum(halo.mass), folder)<1e9:
			ax.add_patch(plt.Circle((center[HALO][0],center[HALO][1]),R200[HALO],fill=False, **kwargs))

def plotDiff(file1, file2, folder,**kwargs):
	
	mmax=7e11
	
	fig=plt.figure()
	ax = fig.add_subplot(1, 1, 1)
	
#	p.diffslice(file1.replace("part","grid"),file2.replace("part","grid"),9,force=0, field="field.d", log=1)

	den, tag, part, star, nstar, a, ngrp = readAll(file1)	
	center   = np.empty((ngrp), dtype=np.object)
	R200     = np.zeros(ngrp, dtype=np.float32)
	for HALO in range(ngrp):
		halo,center[HALO],R200[HALO]   = getHalo(part,tag,den,HALO)	
		if physique.m2mo(np.sum(halo.mass), folder)<mmax:
			ax.add_patch(plt.Circle((center[HALO][0],center[HALO][1]),R200[HALO],fill=False, color='b'))
	plt.plot(star.x, star.y, 'b.', label="wSN")

	den, tag, part, star, nstar, a, ngrp = readAll(file2)	
	center   = np.empty((ngrp), dtype=np.object)
	R200     = np.zeros(ngrp, dtype=np.float32)
	for HALO in range(ngrp):
		halo,center[HALO],R200[HALO]   = getHalo(part,tag,den,HALO)	
		if physique.m2mo(np.sum(halo.mass), folder)<mmax:
			ax.add_patch(plt.Circle((center[HALO][0],center[HALO][1]),R200[HALO],fill=False, color='r',**kwargs))
	plt.plot(star.x,star.y, 'r.', label="woSN")
	
	plt.legend()
	plt.show(block=False)

def plot_bkp(file):
		
	den, tag, part, star, nstar, a, ngrp = readAll(file)	
	#tree = kd.Tree(star.x,star.y,star.z)

	center   = np.empty((ngrp), dtype=np.object)
	R200     = np.zeros(ngrp, dtype=np.float32)
	mask     = np.empty((ngrp), dtype=np.object)

	for HALO in range(ngrp):

		halo,center[HALO],R200[HALO]   = getHalo(part,tag,den,HALO)
		#starHalo, mask[HALO]           = findStar(center[HALO],R200[HALO],star, tree)

		#plt.plot(starHalo.x[HALO],starHalo.y[HALO], '.')
		#plt.plot(center[HALO][0],center[HALO][1], 'ro')	
		plt.Circle((center[HALO][0],center[HALO][1]),R200[HALO],fill=False)

	plt.show(block=False)
	


def massDistri(Mh,folder,**kwargs):
	
	#watson
	Nbins=16
	m=physique.m2mo(Mh,folder)
	
	#m_min=np.min(m)
	#m_max=np.max(m)
	
	m_min=173057153.955
	m_max=12656926065.8

	
	bins = np.linspace(m_min,m_max,num=Nbins+1)
	n,bin_edges=np.histogram(m,bins=bins)
	
	bin_mean=(bin_edges[1:]+bin_edges[:-1])/2		
	
	return np.histogram(m,bins=bins)
	"""
	plt.loglog(bin_mean,n,'o',**kwargs)
	#plt.plot(np.log10(bin_mean),np.log10(n),'o',**kwargs)
	plt.xlabel(r'DM Mass of Halo ($M_{\odot}$)')
	plt.ylabel(r'#')		
	plt.legend()
	plt.show(block=False)
	"""

