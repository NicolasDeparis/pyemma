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
import luminosity as lum
import sfr
import convert 
import amr

def getStarburstNphot():
	data = np.loadtxt('python/Starburst99/912_inst_e.dat.txt',unpack=True)
	x=data[0]
	y=np.power(10,data[1])
	return x,y

class Halo:
	def __init__(self,ngrp):
		self.t=0
		self.file=""
		
		#self.step=
		self.folder=""
		
		self.N=ngrp
		self.x=np.zeros(ngrp, dtype=np.float32)
		self.y=np.zeros(ngrp, dtype=np.float32)
		self.z=np.zeros(ngrp, dtype=np.float32)
		self.R=np.zeros(ngrp, dtype=np.float32)
		self.mass=np.zeros(ngrp, dtype=np.float32)
		self.starmass=np.zeros(ngrp, dtype=np.float32)
		self.part=np.empty(ngrp, dtype=np.object)
		self.star=np.empty(ngrp, dtype=np.object)

	def getN(self):
		return self.N
	def getM(self):
		return self.mass
	def getMstar(self):
		return self.starmass
	def toFile(self,name):
		with open(name, 'wb') as output:		
			pickle.dump(self, output, pickle.HIGHEST_PROTOCOL)	
					
	def getCubeLimit(self,n):
		R=self.R[n]
		xmin=self.x[n] - R
		xmax=self.x[n] + R
		ymin=self.y[n] - R
		ymax=self.y[n] + R
		zmin=self.z[n] - R
		zmax=self.z[n] + R
		return (xmin,xmax,ymin,ymax,zmin,zmax)

	def getfield(self,file,n,field,**kwargs):
		bound = getCubeLimit()
		level=7
		Map=p.getSlice(file,level, xmin=xmin, xmax=xmax ,ymin=ymin, ymax=ymax ,zmin=zmin, zmax=zmax , **kwargs)
		return bound,Map
	def extractCube(self,data,halo_number,level):
		bound= self.getCubeLimit(halo_number)
		xmin,xmax,ymin,ymax,zmin,zmax = bound
		X0=np.linspace(0,1,num=2**level+1)[:-1]	
		x= np.where( (X0>=xmin)&(X0<=xmax) )
		y= np.where( (X0>=ymin)&(X0<=ymax) )
		z= np.where( (X0>=zmin)&(X0<=zmax) )
		X, Y, Z = np.meshgrid(x, y, z)
		
		if (xmin>=0)&(ymin>=0)&(zmin>=0)&(xmax<=1)&(ymax<=1)&(zmax<=1) :
			print "true"
			cube = data[X,Y,Z]
		else :
			print "false"
			cube = None
		#TODO fix this function	
		return cube 
		
	def mean_r_to_center(self,ihalo):
		#xc = self.x[ihalo]
		#yc = self.y[ihalo]
		#zc = self.z[ihalo]
		
		xc = np.mean(self.star[ihalo].x)
		yc = np.mean(self.star[ihalo].y)
		zc = np.mean(self.star[ihalo].z)
		
		n = len(self.star[ihalo].x)
		
		R = np.zeros(n)
		for ipart in range(n):
			x=self.star[ihalo].x[ipart]
			y=self.star[ihalo].y[ipart]
			z=self.star[ihalo].z[ipart]
			
			
			
			dx2= np.power(xc-x,2)
			dy2= np.power(yc-y,2)
			dz2= np.power(zc-z,2)
			
			
			R[ipart] = np.sqrt(dx2+dy2+dz2)
			
		return np.mean(R)
		
	def getDMinR200(self,n):
		x=self.part[n].x
		y=self.part[n].y
		z=self.part[n].z
		tree=kd.Tree(x,y,z)

		center=(self.x[n],self.y[n],self.z[n])
		mask,dump = tree.locatenear(center,self.R[n])
		part_R200 = part.Part(len(mask),1)
		part_R200.mask(self.part[n],mask)
		
		Mdm = np.sum(part_R200.mass)
		return Mdm
		

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
		
	def plot3d(self,x,y,z,*args,**kwargs):
		"""
		Simple 3d plotting function. If a figure is provided through
		the keyword "fig", it is used, and no plt.show() will be
		performed inside the function:
		   plot3d( x, y, z, fig = figure12 )
		"""

		if kwargs.has_key("fig"):
			fig =  kwargs.pop("fig")
			Show = False
		else:
			fig = plt.figure()
			Show = True
		ax=fig.add_subplot(111,projection="3d")
		ax.plot(x,y,z,'.',*args,**kwargs)
		if Show:
			plt.show()
		return ax

	def getNphot(self,n,t0):
		N_dot=0
		t,Nphot = getStarburstNphot()
		
		for age in self.star[n].age:				
			T=t0-age
			photStar=Nphot[(np.abs(t-T)).argmin()]/1e50
			N_dot+=photStar
		return N_dot*1e50

	def getLuminosity(self,n):
		modelmag, modelage = lum.getModel()		
		return lum.sumMag(self.star[n].mass, self.star[n].age, modelmag, modelage)

	def SFR(self,n,folder,**kwargs):
		z,fr = sfr.getFromPop(self.star[n],folder)
		plt.plot(z,fr,**kwargs)
		plt.yscale("log")
		plt.legend()
		
####################################################################################################
####################################################################################################
	
def read(name):	
	try :
		with open(name, 'rb') as input:
			return pickle.load(input)			
	except  IOError:
		genFiles(name)
		
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
	return np.power(3.*M/(200.*4.*np.pi),1./3.)

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
	"""compute kd tree if it doesn't already exist, else read it"""
	name=file.replace("part.","halo.")+ ".%stree"%type
	try :
		open(name)
		tree=kd.Tree.FromFile(name)
	except  IOError:
		tree=kd.Tree(part.x,part.y,part.z)
		tree.WriteToFile(name)	
	return tree

def genFiles(file):
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
			halo.mass[iHalo]=np.sum(part_halo.mass)
			halo.starmass[iHalo]=np.sum(star_halo.mass)

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
	halo = read(name)
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
		
	#print nstar_halo/mtot_star
	plt.title("%s\t%s%s"%(file,axe1,axe2))
	plt.xlim(0,1)
	plt.ylim(0,1)

def plotOne(file,	halo_number,axe1="x",axe2="y", **kwargs):
	
	name=file.replace("part.","halo.")+ ".halo"
	halo = read(name)

	nstar,a,star_all = part.read(file.replace("part","star"))
	tree=genTree(file,star_all,"star")
	
	center=(halo.x[halo_number],halo.y[halo_number],halo.z[halo_number])
	R200=halo.R[halo_number] *1
	print file, "R200=%e"%R200
	
	star_halo = getStar(center,R200,star_all,tree)

	fig=plt.figure()
	plt.title(file)
	ax = fig.add_subplot(1,1,1)


	halo.plot(ax,halo_number,axe1,axe2)
	halo.plotNumber(ax,halo_number,axe1,axe2,)
	
	halo.plotPart(halo_number,axe1,axe2,c='k')
	plotStar(star_halo,axe1,axe2,c='r')
	halo.plotStar(halo_number,axe1,axe2,c='g')


def getfield(step,halo_number,field, **kwargs):
	
	file = step.part_path
	name=file.replace("part.","halo.")+ ".halo"
	halo = read(name)
	
	return halo.getfield(step.grid_path,halo_number,field, force=1)

def barionique_fraction(step, **kwargs):

	name=step.part_path.replace("part.","halo.")+ ".halo"
	halo = read(name)

	level=7
	force=0
	field="field.d"

	filename = step.grid_path
	nproc=IO.getNproc(filename)	
	convert.oct2grid(filename,level, force, nproc, field)
	data = amr.cube(filename.replace("grid.","cube"+ str(level)+ ".") + "." + field)
	data = data.getData()
	
	n=halo.getN()
	Mdm = np.zeros(n)
	Mb = np.zeros(n)
	for halo_number in range(n):
		
		bound= halo.getCubeLimit(halo_number)
		xmin,xmax,ymin,ymax,zmin,zmax = bound
		cube= halo.extractCube(data,halo_number,level)

		dx=xmax-xmin
		dy=ymax-ymin
		dz=zmax-zmin
		v = dx*dy*dz
		Mb[i] = np.sum(cube)*v
		#if v:
		#	print dx,dy,dz, cube
		#	Mb[i] = np.max(cube)
		
		Mdm[i] = halo.getDMinR200(halo_number)

	bf = Mb/Mdm
		
	Mdm = physique.m2mo(Mdm,step.step_path)
	plt.loglog(Mdm,bf,'.', **kwargs)
		
	plt.xlabel(r"Halo Mass $M_{\odot}$")
	plt.ylabel(r'Baryonic fraction')
	
def rho_max(step, **kwargs):

	name=step.part_path.replace("part.","halo.")+ ".halo"
	halo = read(name)

	level=7
	force=0
	field="field.d"

	filename = step.grid_path
	nproc=IO.getNproc(filename)	
	convert.oct2grid(filename,level, force, nproc, field)
	data = amr.cube(filename.replace("grid.","cube"+ str(level)+ ".") + "." + field)
	data = data.getData()
	
	n=halo.getN()
	Mdm = np.zeros(n)
	rho_max = np.zeros(n)
	for halo_number in range(n):
		
		bound= halo.getCubeLimit(halo_number)
		xmin,xmax,ymin,ymax,zmin,zmax = bound
		cube= halo.extractCube(data,halo_number,level	)
	
		if cube!=None:
			rho_max[halo_number] = np.max(cube)

		Mdm[halo_number] = halo.getDMinR200(halo_number)

		
	Mdm = physique.m2mo(Mdm,step.step_path)
	plt.loglog(Mdm,rho_max,'.', **kwargs)
		
	plt.xlabel(r"Halo Mass $M_{\odot}$")
	plt.ylabel(r'Maximal density')
	#TODO profile radial du gas
		
def distance_au_centre(step, **kwargs):

	name=step.part_path.replace("part.","halo.")+ ".halo"
	halo = read(name)

	level=7
	force=0
	field="field.d"

	filename = step.grid_path
	nproc=IO.getNproc(filename)	
	convert.oct2grid(filename,level, force, nproc, field)
	data = amr.cube(filename.replace("grid.","cube"+ str(level)+ ".") + "." + field)
	data = data.getData()
	
	n=halo.getN()
	Mdm = np.zeros(n)
	r = np.zeros(n)
	

	for halo_number in range(n):
		
		r[halo_number]=halo.mean_r_to_center(halo_number)
		Mdm[halo_number] = halo.getDMinR200(halo_number)


	
	r=physique.Cell2Meter(r,step.step_path,7)
	Mdm = physique.m2mo(Mdm,step.step_path)
	plt.loglog(Mdm,r,'.', **kwargs)
		
		
	
	mask = np.where(r>0)
	
	logx=np.log10(Mdm[mask])
	logz=np.log10(r[mask])
	b=10
	w,bin_edges=np.histogram(logx, bins=b)
	sfr,bin_edges=np.histogram(logx, bins=b, weights = logz)
	sfr/=w
	sfr=np.power(10,sfr)
	
	bins=(bin_edges[1:]+bin_edges[:-1])/2
	bins=np.power(10,bins)
	
	plt.loglog(bins, sfr)
		
		
	plt.xlabel(r"Halo Mass $M_{\odot}$")
	plt.ylabel(r'Average istance to center $Pc$')

	#TODO profile radial du gas
		
def compare_field(step1,step2,halo_number,axe1="x",axe2="y", **kwargs):
	
	field = "field.d"
	
	bound,f1 = getfield(step1,halo_number,field, force=1)
	bound,f2 = getfield(step2,halo_number,field, force=1)
	
	
	#plt.imshow(np.log10(f2-f1), interpolation='nearest')
	plt.imshow(f2-f1, interpolation='nearest')

	plt.colorbar()

	plt.show()

def plotOne3D(file,	halo_number,axe1="x",axe2="y", **kwargs):
	name=file.replace("part.","halo.")+ ".halo"
	halo = read(name)

	nstar,a,star_all = part.read(file.replace("part","star"))
	tree=genTree(file,star_all,"star")
	
	center=(halo.x[halo_number],halo.y[halo_number],halo.z[halo_number])
	R200=halo.R[halo_number]
	
	star_halo = getStar(center,R200,star_all,tree)

	fig1=plt.figure()
	#ax = fig.add_subplot(1,1,1)

	"""
	x=halo.part[halo_number].x
	y=halo.part[halo_number].y
	z=halo.part[halo_number].z
	halo.plot3d(x,y,z,fig=fig1,color='r')
	"""
	
	x=star_halo.x
	y=star_halo.y
	z=star_halo.z
	halo.plot3d(x,y,z,fig=fig1,color='g')
		


	plt.show()


def compare(file1,file2,halo_number, **kwargs):
	
	name=file1.replace("part.","halo.")+ ".halo"
	halo1 = read(name)
	
	name=file2.replace("part.","halo.")+ ".halo"
	halo2 = read(name)
	
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
	halo = read(name)
	
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
	
	halo = read(file.replace("part.","halo.")+ ".halo")
	n=halo.getN()
	
	mh=np.zeros(n, dtype=np.float64)
	lum=np.zeros(n, dtype=np.float64)
	
	for iHalo in range(n):
		par=halo.part[iHalo]
		mh[iHalo]=np.sum(par.mass)*u_mass
		lum[iHalo] = halo.getLuminosity(iHalo)
	
	size = p["box_size_Mpc/h"]
	V = size**3
	#plt.loglog(mh,lum/V,'.')

	
	Nbins=len(mh)
	bins_edges = np.linspace(np.min(mh),np.max(mh),Nbins+1)
	
	lum,bin_edges=np.histogram(mh,bins=bins_edges)
	
	lum/V/np.diff(bin_edges)
	mh=(bin_edges[1:]+bin_edges[:-1])/2

	
	print mh,lum
	plt.loglog(mh,lum,'o')
	
	
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

def findCenterPop(x,y,z):
	"""
		return the center of mass of a population
	"""
	xc = np.mean(x) 
	yc = np.mean(y)
	zc = np.mean(z)
	return xc,yc,zc
	

def identify(file1, file2, n):
	"""
		locate particle of halo in another snap
	"""
	
	halo1 = read(file1)
	part1 = halo1.part[n]

	N,a,part2 = part.read(file2)
	
	idx=np.int32(part1.idx)
	
	sorted_idx = np.argsort(part2.idx)
	
	x= part2.x[sorted_idx[idx]]
	y= part2.y[sorted_idx[idx]]	
	z= part2.z[sorted_idx[idx]]
	
	center = findCenterPop(x,y,z)
	
	return center
	
def follow(file1, file2, n, **kwargs):

	halo1 = read(file1.part_path)
	halo2 = read(file2.part_path)
	
	R = halo1.R[n]
	
	a1 = part.geta(file1.part_path+".p00000")
	a2 = part.geta(file2.part_path+".p00000")
	
	c1 = findCenterPop(halo1.part[n].x, halo1.part[n].y, halo1.part[n].z)	
	c2 = identify(file1.part_path, file2.part_path, n)

	print a1, c1
	print a2, c2
	
	return R, a1, a2, c1, c2
	
	"""
	s=5.
	
	
	fig=plt.figure()
	ax = fig.add_subplot(1,1,1)
	print xmin, xmax, ymin, ymax, zmin, zmax
	p.slice(file1.grid_path,11, field="rfield.temp", xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax, log=True, force= True)



	
	fig=plt.figure()
	ax = fig.add_subplot(1,1,1)
	print xmin, xmax, ymin, ymax, zmin, zmax
	p.slice(file2.grid_path,11, field="rfield.temp", xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax, log=True, force= True)
	"""

def box(c,R):
	xmin=c[0]-R
	xmax=c[0]+R
	ymin=c[1]-R
	ymax=c[1]+R
	zmin=c[2]-R
	zmax=c[2]+R
	return xmin, xmax, ymin, ymax, zmin, zmax
	
def interpPos(a, a1, a2, c1, c2):
	x= np.interp(a, (a1,a2), (c1[0], c2[0]))
	y= np.interp(a, (a1,a2), (c1[1], c2[1]))
	z= np.interp(a, (a1,a2), (c1[2], c2[2]))
	return x, y, z
	
def followRun(folder, file1, file2, n, **kwargs):
	
	R, a1, a2, c1, c2 = follow(file1, file2, n, **kwargs)
	

	listdir = os.listdir(folder)
	mask = np.where
	
	a = (a1+a2)/2	
	interp_c = interpPos(a, a1, a2, c1, c2)
	
	xmin, xmax, ymin, ymax, zmin, zmax = box(interp_c, R*5)

	fig=plt.figure()
	p.slice(file1.grid_path,11, field="rfield.temp", xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax, log=True, force= True)
