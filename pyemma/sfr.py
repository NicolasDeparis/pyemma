import numpy as np
import matplotlib.pylab as plt

import physique
import part
import IO

def fromSnap(file, lab):

	folder,filename = IO.splitPath(file)
	folder +="/"

	param = IO.Params(folder = folder).get()
	L = float(param["unit_l"])/3.085677e16 /1e6
	print L

	Ntot,a,stars=part.read(file)
	
	if len(stars.mass):
		b=8
		n0,bin_edges=np.histogram(stars.age,bins=b)
	
		z=physique.a2z(physique.t2a( bin_edges))

		n=len(n0)
		sfr=np.zeros(n)
		for i in range(n):
			sfr[i] = physique.m2mo( stars.mass[0] * (n0[i]) ,folder )  / float( bin_edges[i+1] - bin_edges[i]) / pow(L,3)
	
		plt.semilogy(z[:-1],sfr, label = lab)
		
		plt.xlabel('z')
		plt.ylabel(r'$M_{0}.yrs^{-1}.Mpc^{-3}.h$')	
		plt.legend()
		plt.show(block=False)

	
def observation():
	
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
