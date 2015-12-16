import os, sys
import numpy as np
import matplotlib.pylab as plt

import amr
import convert
import step
import sedov

from math import sqrt

def compprofile(data=None):
	
	
	"""
	#Comparaison resolution
	l6=step.step(9,"data/thcell_l6/data/")
	l7=step.step(9,"data/thcell_l7/data/")
	l8=step.step(10,"../data_hpc/thcell_8/data/")
	
	plotProfile(l6,6,"L6 ", force=0)
	plotProfile(l7,7,"L7 ", force=1)
	plotProfile(l8,8,"L8 ", force=0)
	"""
	
	"""
	#comparaison raffinement
	l6=step.step(10,"data/thcell_l6/data/")
	l6p1=step.step(10,"data/thcell_l6p1/data/")
	l7=step.step(10,"data/thcell_l7/data/")
	
	plotProfile(l6,6,"L6 ", force=0)
	plotProfile(l7,7,"L7 ", force=0)
	plotProfile(l6p1,7,"L6+1 ", force=0)
	"""
	
	"""
	#Comparaison raffinement Zoom
	l6p2=step.step(1,"data/thcell_l6p2/data/")
	l6p2_zoom=step.step(1,"data/thcell_l6p2_zoom/data/")
	l8=step.step(2,"../data_hpc/thcell_8/data/")	
	
	plotProfile(l6p2,8,"L6p2 ", force=0)
	plotProfile(l6p2_zoom,8,"L6p2_zoom ", force=0)
	plotProfile(l8,8,"L8 ", force=0)
	"""
	
	#Comparaison thermique cinetique
	#the=step.step(7,"../data_hpc/thoct_l6p4_zoom/data/")
	#the=step.step(7,"../data_hpc/thcell_l6p4_zoom/data/")
	#the=step.step(7,"data/thcell_l6p4_zoom/data/")
	#kin=step.step(6,"data/kin_l6p4_zoom/data/")
	#the=step.step(7,"../data_hpc/kin_l6p4_zoom/data/",injection_type="oct")
	#the=step.step(12,"data/",injection_type="cell")
	#the=step.step(10,"data/thcell_8/data/",injection_type="cell")
	



	#plotProfile(kin,10,"Kinetic ", force=0)	
	
	#axes = ["radial"]
	axes = ["axial","diag_cube","diag_slice","radial" ]
	#axes = ["axial","diag_cube","diag_slice" ]
	
	
	
	#fields = ["field.d", "field.p", "field.vel"]
	field="field.vel"
	i=0
	plt.figure(1,figsize=(15,15))

	data = loaddata(the,10, field=field)

	sedov.gen_solution(the)
	for axe in axes :	
	#for field in fields:
		plt.subplot(2,2,i)
		plt.title(axe)

		plotProfile(the,10,"Thermal ", data=data,force=0,log=1, axe=axe,field=field)
		#plt.ylabel('overdensity')	
		sedov.plot(field)
		i+=1
		
def compprofile_field(data=None):
	
	profiles=[]
	
	#profiles.append([step.step(1,"../data_hpc/thcell_l6p4_zoom/data/",injection_type="cell"), "Thermal cell"])
	#profiles.append([step.step(4,"../data_hpc/thcell_l6p4_zoom/data/",injection_type="cell"), "Thermal cell"])
	#profiles.append([step.step(7,"../data_hpc/thcell_l6p4_zoom/data/",injection_type="cell"), "Thermal cell"])

	#profiles.append([step.step(7,"../data_hpc/thoct_l6p4_zoom/data/",injection_type="oct"),"Thermal oct"])
	#profiles.append([step.step(7,"../data_hpc/kin_l6p4_zoom/data/",injection_type="oct"),"Kinetic"])
	
	#thcell=step.step(1,"../data_hpc/comparaison_ramses/data/",injection_type="cell")
	#thcell=step.step(1,"../data_hpc/thcell_l6p4_zoom/data/",injection_type="cell")
#	thcell=step.step(3,"../data_hpc/comparaison_ramses/data/",injection_type="cell")
	
#	profiles.append([step.step(10,"../data_hpc/thcell_zoom_1Mpc/data/",injection_type="cell"),""])
	#profiles.append([step.step(1,"../krusty/Quartz/data/test/data/",injection_type="cell"),""])
	#profiles.append([step.step(4,"../krusty/Quartz/data/test2/data/",injection_type="cell"),""])
	
	
	profiles.append([step.step(8,"../krusty/Quartz/data/test3/data/",injection_type="cell"),""])
	



	
	axe = "radial"

	#axes = ["axial","diag_cube","diag_slice","radial" ]
	#axes = ["axial","diag_cube","diag_slice" ]
		
	fields = [	["field.d","Density"],
				["field.p","Pressure"],
				["field.vel","Velocity"] ]
	
	#field="field.d"
#	plt.figure(1,figsize=(15,5))
	
	

	
	
	label="Thermal cell"
	#sedov.gen_solution(profiles[0][0])
	trig=0
	for profile in profiles:
		i=0		
		sedov.gen_solution(profile[0])
		
		for field in fields:
			j=0
			#plt.subplot(1,3,i)
			plt.figure(i)
			plt.title(field[1])
		
			#if trig == 0:
			if True:
				sedov.plot(field[0])
			
			plotProfile(profile[0],10,profile[1], data=data,force=0,log=0, axe=axe,field=field[0])

			#if trig == 0:
			if True:
				plt.legend()
		
			i+=1
			j+=1
		trig=1
		plt.legend()
		
def plotProfile(step,level,lab,log=0,data=None,axe="radial",field="field.u",error=0,**kwarg):
	
	if axe == "axial":
		#axis=["x","y","z"]
		axis=["x"]
		for axe in axis:		
			Nx,Ny,Nz,t,x,y = getProfile(step,level,data=data,axis=axe,**kwarg)
			
	if axe == "diag_cube":
		Nx,Ny,Nz,t,x,y = getProfileDiag(step,level,data=data,**kwarg)
	if axe == "diag_slice":
		Nx,Ny,Nz,t,x,y = getProfileDiag(step,level,data=data,type=1,**kwarg)
	if axe == "radial":
		t,x,y,xerr,yerr = getProfileProper_alloct(step,level,field=field,	**kwarg)
			
	
		
	fig, ax1 = plt.subplots()
	
	error=1
	if error:
		ax1.errorbar(x,y, c='r',xerr=xerr, yerr=yerr, label="Numerical")
	else:
		if log:
			ax1.semilogy(x,y, label="Numerical")
		else:
			ax1.plot(x,y,'o-', label="Numerical")

	"""
	data =  amr.allOct(step.grid_path,field=field,**kwarg)
	center = get_center("center",step.injection_type,level)
	R=getR_alloct(data, center)	
	rhos=data.getMap()
	ax1.scatter(R,rhos)
	"""
			
	
	sedov.plot(ax1, field)
	ax1.legend()


	ax2 = ax1.twinx()
	t,x,y,xerr,yerr = getProfileProper_alloct(step,level,field="level",	**kwarg)
	ax2.plot(x,y, "--k",label="Level")
	ax2.set_ylabel('level')
	ax2.set_ylim(6,11)
	ax2.legend(loc=4)


	

	#plt.ylim(0,3)
	#plt.xlabel('x [box unit]')
	


def plot(step,field, title, axe, level, force=0, log=0):
    t=sedov.gen_solution(step[0])
    plotProfile(step[0],level,step[1], force=force,log=log, axe=axe,field=field)    
    #plt.title(title + " at t=%e"%t)
    #plt.legend()
    

def loaddata(step,level,field="field.d",**kwarg):	
	fileName = step.grid_name
	convert.oct2grid(step.grid_path,level,field=field,**kwarg)
	file = amr.cubename(fileName,level,field)
	data =  amr.cube(step.grid_folder+file)
	return data


def getProfile(step,level, data=None, axis="x",**kwarg):
	
	if data==None:
		data = loaddata(step,level)

	Nx,Ny,Nz = data.getSize()
	c = data.getData()
	t = data.geta()
	print t
	
	x = np.arange(Nx/2, dtype=np.float32)

	if axis == "x":
		y = c[Nx/2:,Ny/2,Nz/2]
	if axis == "y":
		y = c[Nx/2,Ny/2:,Nz/2]
	if axis == "z":
		y = c[Nx/2,Ny/2,Nz/2:]
		
#	x = np.arange(Nx, dtype=np.float32)
#	y = c[0:Nx,0,0]
	
	x/=Nx
	
	return Nx,Ny,Nz,t,x,y

def getProfileDiag(step,level,data=None, field="field.d", axis="x", type=0,**kwarg):
	
	if data==None:
		data = loaddata(step,level)
		
	c = data.getData()
	t = data.geta()

	Nx,Ny,Nz = data.getSize()

#	x = np.arange(Nx/2)
#	y = c[Nx/2:Nx,Ny/2,Nz/2]

	center = 1
	
	if center:		
		
		if type ==0:
			#diag slice 
			x = np.arange(Nx/2, dtype=np.float32) * np.sqrt(2)
			y = np.zeros(Nx/2)
			for i in range(Nx/2):
				y[i] = c[Nx/2+i,Ny/2+i,Nz/2]

		if type ==1:
			#diag cube
			x = np.arange(Nx/2, dtype=np.float32) * np.sqrt(3)
			y = np.zeros(Nx/2)
			for i in range(Nx/2):
				y[i] = c[Nx/2+i,Ny/2+i,Nz/2+i]
	else:
		x = np.arange(Nx, dtype=np.float32) * np.sqrt(3)
		y = np.zeros(Nx)
		for i in range(Nx):
			y[i] = c[i,i,i]
				
	x/=Nx

	return Nx,Ny,Nz,t,x,y

def getProfileProper(step,level,data=None, xc=0,yc=0,zc=0, R=None,field="field.d", **kwarg):
	
	if data==None:
		data = loaddata(step,level)
	
	Nx,Ny,Nz = data.getSize()
	xc,yc,zc = Nx/2,Ny/2,Nz/2 	 	
 	
	if R==None:
		print "getting r"
		
		x=np.arange(Nx)
		y=np.arange(Ny)
		z=np.arange(Nz)
		
		X=x-xc; x2 = np.power(X,2)
		Y=y-yc; y2 = np.power(Y,2)
		Z=z-zc; z2 = np.power(Z,2)
		
		grid = np.meshgrid(x, y, z, indexing='ij')
		r2 = x2[grid[0]]+y2[grid[1]]+z2[grid[2]]
		R=np.int32(np.sqrt(r2))
	

	Nbins = int(Nx+Ny+Nz)/3	
	Nbins /= 4
	print "nbins=",Nbins
	
	MIN=np.min(R)
	MAX=np.max(R)
	bin_edges = np.linspace(MIN,MAX,Nbins+1)
	
	
	w=data.getData()
	
	n,_= np.histogram(R, bins=bin_edges)
	rho,_= np.histogram(R, bins=bin_edges, weights=w)	
	rho2,_= np.histogram(R, bins=bin_edges, weights=w*w)	
	
	
	x=(bin_edges[1:]+bin_edges[:-1])/2
	y = rho/n
	yerr = np.sqrt(rho2/n - y*y)  / np.sqrt(n)*3   
	
	
	t = data.geta()
	return Nx,Ny,Nz,t,x,y, yerr


def get_center(source_pos,feedback_type, level):

	if source_pos=="center":
		xc=0.5
		yc=0.5
		zc=0.5
	if source_pos=="0":
		xc=0.
		yc=0.
		zc=0.
	
	
	dx = np.power(0.5,level)
	
	if feedback_type=="cell":
		xc+=dx/2.
		yc+=dx/2.
		zc+=dx/2.
	if feedback_type=="oct":				
		xc+=dx
		yc+=dx
		zc+=dx
	
	
	return xc,yc,zc

def getR_alloct(data, center):
	
	dxs=np.power(0.5,data.getLevel()+1)
	#dxs *=0
	
	x= np.add(data.getx(),dxs)
	y= np.add(data.gety(),dxs)
	z= np.add(data.getz(),dxs)

	x2=np.power(np.subtract(x,center[0]),2.)
	y2=np.power(np.subtract(y,center[1]),2.)
	z2=np.power(np.subtract(z,center[2]),2.)
	
	return np.sqrt(x2+y2+z2)
	
def getProfileProper_alloct(step,level,data=None,field="field.d", src_pos="center", **kwarg):

	if data==None:
		data =  amr.allOct(step.grid_path,field=field,**kwarg)

		
	center = get_center(src_pos,step.injection_type,7)
	R=getR_alloct(data, center)
	
	"""
	if src_pos=="center":
		R =R/2.
	"""
	
	
	rhos=data.getMap()
	
	#print "posmax",R[np.argmax(rhos)]


	"""
	mask=np.where(rhos>1)
	plt.plot(R[mask],rhos[mask],'.', label="all points")
	"""

	"""
	plotSedov(step,level, t=data.geta())
	plt.xlim(0,0.15)
	"""
	
	nbins = 2**(level-1)
	"""
	if src_pos == "center":
		nbins = 2**(level-1)
	if src_pos == "0":
		nbins = 2**(level)
	"""
	
	rho2,   bins = np.histogram(R, bins=nbins, weights=rhos*rhos)
	rho1,   bins = np.histogram(R, bins=nbins, weights=rhos)
	n,  	bins = np.histogram(R, bins=nbins)

	"""
	dr = 0.02
	rc = 0.01
	mask = np.where((R>rc-dr) & (R<rc+dr) & (data.getz()>0.499) & (data.getz()<0.501))
	
	plt.figure()
	plt.scatter(data.getx()[mask],data.gety()[mask],c=np.log(R[mask]), edgecolors='None')
	#plt.plot(data.getx()[mask],data.gety()[mask], '.' )
	plt.plot(center[0],center[1],'o')
	exit
	"""
	
	xplot=(bins[1:]+bins[:-1])/2.
	yplot = np.divide(rho1,n)


#n, _ = np.histogram(x, bins=nbins)
#sy, _ = np.histogram(x, bins=nbins, weights=y)
#sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
	
	

	xerr = np.diff(bins)/2.
	yerr = np.sqrt(rho2/n - yplot*yplot) / np.sqrt(n) *3.
	#	plt.plot(xplot,yplot, 'o-', label = "numerical")

	return data.geta(),xplot,yplot, xerr, yerr
		
	
def getR(Nx,Ny,Nz, xc,yc,zc):
	print "getting r"
	
	x=np.arange(Nx)
	y=np.arange(Ny)
	z=np.arange(Nz)
	
	X=x-xc; x2 = np.power(X,2)
	Y=y-yc; y2 = np.power(Y,2)
	Z=z-zc; z2 = np.power(Z,2)
	
	grid = np.meshgrid(x, y, z, indexing='ij')
	r2 = x2[grid[0]]+y2[grid[1]]+z2[grid[2]]
	r=np.int32(np.sqrt(r2))
	
	return r

	

	
def comparField(args):
	args.field = ["field.d"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="density")

	args.field = ["field.p"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="P")

	args.field = ["field.u"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="Vx")
	
def comparMethod(args):
	args.field = ["field.d"]
	Nx,Ny,Nz,x,y = getProfile(fileName,args)
	plt.plot(x,y, label="x")
	
	Nx,Ny,Nz,x,y = getProfileDiag(fileName,args)
	plt.plot(x,y, label="diag")
	
	Nx,Ny,Nz,x,y = getProfileProper(fileName,args)
	plt.plot(x,y, label="proper")

def getAllProfil(args, folder = "data/", force=0, field='field.d'):

	Nx,Ny,Nz = 128,128,128
	xc,yc,zc = 64,64,64
	r = getR(Nx,Ny,Nz,xc,yc,zc)

	files = np.sort(os.listdir(folder))
	for file in files:
		if "grid." in file and ".p00000" in file: 
			if file.replace("grid","profile")[:-7] in files:
				print file.replace("grid","profile") + " allready exist"
				continue
			file= "%s%s"%(folder,file[:-7])

			convert.oct2grid(file,args.level,force=force, field=field)
			Nx,Ny,Nz,x,y = getProfileProper(amr.cubename(file,args.level,field), R=r)				
			np.savetxt(file.replace("grid","profile"),(x,y))			

def readAllProfil(folder="data/"):
	data = []
	for file in np.sort(os.listdir(folder)):		
		if "profile." in file : 		
			data.append( np.loadtxt("%s%s"%(folder,file)) )
	return np.array(data)

def plotAllProfil(data):
	for set in data:
		plt.plot(set[0],set[1])
	plt.show()

def findFrontPosition(data):
	position = np.zeros(len(data))
	i=0
	for set in data:
		grad = np.gradient(set[1])

		xinterp = np.linspace(0,len(data),10000)
		yinterp = np.interp(xinterp, set[0], grad)

		signchange = (np.diff(np.sign(yinterp))!=0)*1
		wheresignchange = np.where(signchange==1)
		
		if np.size(wheresignchange):			
			ind = np.argmax(set[1][np.int_(xinterp[wheresignchange[0]])]) #peut etre revoir l'arrondis ??
			pos = wheresignchange[0][ind]
			#pos = np.int_(xinterp[wheresignchange[0][ind]])
			#pos = xinterp[ind]
				
			#pos = set[0][np.int_(xinterp[wheresignchange[0][0]])]
			print pos
		else :
			pos = 0
		
		position[i] = xinterp[pos]
		
		i+=1

	return position

def findFrontPositionXion(data):
	position = np.zeros(len(data))
	i=0
	for set in data:
		
		x = set[0]
		y = set[1]
		
		xinterp = np.linspace(0,len(data),4196)
		yinterp = np.interp(xinterp, set[0], y)
		
		signchange = (np.diff(np.sign(yinterp-0.5))!=0)*1
			
		wheresignchange = np.where(signchange==1)
		
		if np.size(wheresignchange):			
			ind = np.argmax(set[1][np.int_(xinterp[wheresignchange[0]])]) #peut etre revoir l'arrondis ??
			pos = wheresignchange[0][ind]
		
			print pos
		else :
			pos = 0
		print wheresignchange

		position[i] = xinterp[pos]
		
		i+=1

	return position

	
def findFrontAmp(data):
	amp = np.zeros(len(data))
	i=0
	for set in data:
		m = [0,0]
		for j in range(2):
			m[j] = np.nanmax(set[1])
			set[1][np.nanargmax(set[1])] = 0
			
		amp[i] = np.nanmax(set[1])
		i+=1
	return amp

def frontSpeed(x,t):	
	return np.gradient(x[1:],np.diff(t))
