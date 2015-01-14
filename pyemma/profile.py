import os, sys
import numpy as np

import amr
import convert

from math import sqrt

def getProfile(fileName,level,field):
	file = amr.cubename(fileName,level,field)
	data =  amr.cube(file)
	Nx,Ny,Nz = data.getSize()	
	c = data.getData()
	
	x = np.arange(Nx/2)
	y = c[Nx/2:Nx,Ny/2,Nz/2]
	
	return Nx,Ny,Nz,x,y

def getProfileDiag(fileName,args):
	level = args.level
	field = args.field[0]

	file  = fileName.replace("grid","cube"+ str(level)) + "." + field

	data =  amr.cube(file)

	c =  data.getData()
	Nx,Ny,Nz = data.getSize()

#	x = np.arange(Nx/2)
#	y = c[Nx/2:Nx,Ny/2,Nz/2]

	center = 0
	
	if center:
		diag = int(np.sqrt(Nx*Nx/4 + Ny*Ny/4))

		x = np.arange(Nx/2) * np.sqrt(2)
		y = np.zeros(Nx/2)
		for i in range(Nx/2):
			y[i] = c[Nx/2+i,Ny/2+i,Nz/2]
	else:
		diag = int(np.sqrt(Nx*Nx + Ny*Ny))

		x = np.arange(Nx) * np.sqrt(2)
		y = np.zeros(Nx)
		for i in range(Nx):
			y[i] = c[i,i,i]
		

	return Nx,Ny,Nz,x,y

def getR(Nx,Ny,Nz, xc,yc,zc):
	print "getting r"
	r = np.zeros(Nx*Ny*Nz)
	r = r.reshape(Nx,Ny,Nz)
	
	for z in range(0,Nz):
		Z = float(z-zc)
		for y in range(0,Ny):
			Y = float(y-yc)
			for x in range(0,Nx):
				X = float(x-xc)

				r[z][y][x] = int(np.sqrt( X*X + Y*Y + Z*Z ) )
	return r

def getProfileProper(fileName, xc=0,yc=0,zc=0, R=None):
	print "getting profile of %s"%fileName
	data = amr.cube(fileName)

	Nx,Ny,Nz = data.getSize()
	xc,yc,zc = Nx/2.,Ny/2.,Nz/2.

	diag = int((Nx-xc) * sqrt(3.0) ) +1
		
	c =  data.getData()
	rho = np.zeros(diag, dtype=np.float)
	n = np.zeros(diag, dtype=np.float)
	 
	if R==None:
		R =getR(Nx,Ny,Nz,xc,yc,zc)

	for z in range(0,Nz):		
		for y in range(0,Ny):			
			for x in range(0,Nx):
				r = R[z][y][x]
				rho[r] += c[z][y][x]
				n[r] += 1

	yplot = np.divide(rho,n)
	xplot = np.arange(diag)

	return Nx,Ny,Nz,xplot,yplot
	
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

	Nx,Ny,Nz = 64,64,64
	xc,yc,zc = 32,32,32
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
