#!/usr/bin/env python
import os, sys
import numpy as np
from pyemma  import *
import matplotlib.pylab as plt
import sedov


def plotAnimProfil(args):
	
	field = "field.d"
	
	args.folder="../krusty/Quartz/data/"
	data = profile.readAllProfil(folder = "%sprofile/%s/"%(args.folder, field))	
	t = amr.gettmap(folder = args.folder)
	unit = physique.Cell2Meter(args)/1000

	args.folder="../data/strom_2e6_SN/" 
	data2 = profile.readAllProfil(folder = "%sprofile/%s/"%(args.folder, field))	
	t2 = amr.gettmap(folder = args.folder)
	unit2 = physique.Cell2Meter(args)/1000

	
	i=0
	folder_out = 'img/'
	for i in range(len(data)):
		plt.clf()
			
		plt.semilogy(np.multiply(data[i][0], unit), data[i][1] , label="with SN")
		plt.semilogy(np.multiply(data2[i][0], unit2),data2[i][1], label="without SN")
		
		plt.ylabel(r'delta rho')
		plt.xlabel(r'R (kpc)')
		plt.title("%s Myr"%"{:.2f}".format(t[i])) 

		plt.ylim(1e-1,3)
		plt.legend()
		plt.savefig(folder_out +str(i))
	
if __name__ == "__main__":	

	args = IO.getargs()
	
#	profile.getAllProfil(args, folder=args.folder, force=1, field='field.d')
	plotAnimProfil(args)

#	y = profile.findFrontPositionXion(data)
#	y = np.multiply(y,1./64)
	
	"""
	plt.plot(t[1:len(y)],y[1:])

	plt.ylabel(r' Rfront/Lbox')
	plt.xlabel(r'time (Myr)')
	
	plt.show()
	"""

