import numpy as np
import IO
import param
import matplotlib.pylab as plt

from scipy import integrate
from scipy import stats

class Constantes : 
	def __init__(self):

		self.Parsec = 3.08567758e16				# parsec in m
		self.G      = 6.67384e-11            	# newton constant 
		self.M0     = 1.9891e30              	# Solar masse
		self.c      = 299792458              	# velocity of light in m/s
		self.proton_mass =1.67262158e-27		#kg
		self.planck =  6.62617e-34				#J/s Planck constant
		self.boltzmann = 1.38066e-23 			#J/K Boltzmann constant
		self.eV = 1.602176565e-19 				#electron Volt in Joules
		
		self.Tyr = 977.8 

		self.H0 = 67.0			     			# Hubble constant 
		self.H0_SI = self.H0*1e3/1e6/self.Parsec # Hubble constant in SI unit
		self.h  = self.H0/100.
	
		self.WB = 0.049
		self.WM = 0.3175                 		# Omega(matter)
		self.WV = 0.6825               			# Omega(vacuum) or lambda
		self.WR = 4.165E-5/(self.h*self.h)    	# Omega(radiation)
		self.WK = 1-self.WM-self.WR-self.WV		# Omega curvaturve = 1-Omega(total)

		self.yr = 31556926						# secondes in 1 year


def rho_c():
	""" critical density """
	c = Constantes()
	return 3.*c.H0_SI**2/(8.*np.pi*c.G)

########################################################################
# time
########################################################################

def a2z(a) :
	return 1.0/a -1.0

def z2a(z):
	return 1.0/(1.0+z)

def a2t(a) :
	"""
		convert expansion factor to time
	"""
	c = Constantes()

	az = 1.0/(1.0+1.0*a2z(a))
	age = 0.
	n=1000         # number of points in integrals
	for i in range(n):
		a = az*(i+0.5)/n
		adot = np.sqrt(c.WK+(c.WM / a)+(c.WR/ (a*a) )+ (c.WV*a*a) )
		age = age + 1./adot
	zage = az*age/n
	zage_Gyr = (c.Tyr/c.H0)*zage

	return zage_Gyr*1e9

def t2a(t):
	"""
		convert time to expansion factor
	"""
	n = 10000
	A = np.arange(n+1) / float(n)
	T = a2t(A)
	return np.interp(t,T,A)

########################################################################
# mass
########################################################################

def m2mo(m,P) :
	"""
		convert code mass in solar mass
	"""
	c = Constantes()	
	p = P.info.get()
	unit_m = float(p["unit_mass"])
	return m/c.M0*unit_m 

def mo2m(mo,P) :
	"""
		convert solar mass in code unit
	"""
	c = Constantes()
	p = P.info.get()
	unit_m = float(p["unit_mass"])
	return mo*c.M0/unit_m 

########################################################################
# length
########################################################################

def Cell2Meter(x,P,level):
	"""
		return the size of a cell in parsec
	"""
	p = P.info.get()
	L = float(p["unit_l"])/3.085677e16
	dx = pow(2,- level)*L
	return x*dx
	

########################################################################
# energy
########################################################################

def E2lambda(E):
	cst=Constantes() 
	c=cst.c
	h=cst.planck
	return h*c/E

def lambda2E(l):
	cst=Constantes() 
	c=cst.c
	h=cst.planck
	return h*c/l
	
def eV2lambda(eV):
	cst=Constantes()
	h = cst.planck
	c = cst.c	  
	Ej = eV*cst.eV
	return h*c/Ej
	
def lambda2nu(l):
	cst=Constantes() 
	return cst.c/l
	
def nu2E(nu):
	cst=Constantes()
	return cst.planck*nu
	

########################################################################
# photon
########################################################################
	
def Nphotons_1():
	#iliev 2006 http://adsabs.harvard.edu/abs/2006MNRAS.369.1625I
	c = Constantes()

	o_b = c.WB
	o_0 = c.WM+c.WR+c.WV
	mp = c.proton_mass
	
	f_star=50000
	f_esc=0.2
	Ni=0.2
	f_gamma = f_star * f_esc *Ni
	
	M=1.
	dti=1.
	
	N_dot = f_gamma*M*o_b/(dti*o_0*mp)
	
	print N_dot


def Nphotons_2():
	#E_rayonne = 100 * E_SN
	SN_EGY = 8.73e11/1.736e-2
	N_SNII = 1.#1.736e-2 
	M = 30 *2e30
	
	E_SN = SN_EGY*N_SNII*M

	Ephot =29.609895722*1.6022e-19
	t_life = 2e7 *31556926
	
	# E_rayonne = phot/s *tlife * Ephot =100 * E_SN
	# phot/s = (100 * E_SN)/ (tlife * Ephot)
	
	print (100*E_SN)/(t_life*Ephot) /M
	
def Nphotons_3():
		#baek 2010 http://adsabs.harvard.edu/abs/2010A%26A...523A...4B
	E = 2.14e45 #erg/s	
	E *= 1e-7 #J/s
	
	Ephot =29.609895722*1.6022e-19
	
	M=np.linspace(8,120,1e5)
	imf=IMF_salpeter(M)
	mean=np.average(M,weights=imf)
	
	print mean
	
	print E/Ephot/(mean*2e30)
	
def Nphotons_4():
	#starburst99 fig77 
	E0 = 10**(52.71)
	c = Constantes()
	M0 =c .M0
	print E0 /(1e6*M0)
	
########################################################################
		
def test():
	M=np.linspace(1,120,1e2)
	plt.loglog(M,BB_int(M))
	

def  tlife_m(M):
	t0=1e10
	return t0*np.power(M,-2.5)

def T_m(M):
	"""temperature function of mass"""
	T0=5800
	alpha=2.5/4
	return T0*np.power(M,alpha)

def L_m(M):
	"""Luminosite function of mass"""
	L0=3.846e26
	alpha=3.5
	return L0*np.power(M,alpha)

def IMF_salpeter(M):
	"""Salpeter Initial Mass Function"""
	E0=1.
	alpha=-2.35
	return E0*np.power(M,alpha)
	
def percent_mass_imf():
	"""compute percent of mass in [6,100]M0 for a salpeter IMF"""
		
	mtot=np.linspace(1,100, 256)
	ntot=IMF_salpeter(mtot)

	m=np.linspace(6,100, 256)
	n=IMF_salpeter(m)
	
	def integ(n,m):
		return integrate.trapz(n*m,m)
	
	print integ(n,m)/integ(ntot,mtot)
	

def black_body(l,T):
	"""Planck law"""
	cst=Constantes()
	h = cst.planck
	c = cst.c
	k = cst.boltzmann
		
	c2 = np.power(c,2)
	l5 = np.power(l,5)
	
	return 2*h*c2/l5   *  1./(np.exp(h*c/(k*l*T))-1.)
	
def BB_int(T):
	lmin=0.
	lmax=np.inf
	return integrate.quad(black_body, lmin, lmax,args=(T,))[0]
	
def Enu(l,M):	
	return L_m(M)* black_body(l,T_m(M)) / BB_int(T_m(M))

def Etot(l,t):
	d=0
	
	
def farUVbelow912():
	data = np.loadtxt('python/Starburst99/912_inst_e.dat.txt',unpack=True)
	x=np.log10(data[0])
	y=data[1]
	plt.plot(x,y)
	
	mask = np.where(x<6)
	A= stats.linregress(x[mask],y[mask])
	x_in=np.linspace(3,9,100)
	y_in=A[0]*x_in+A[1]
	print A
	plt.plot(x_in,y_in)
		
	mask = np.where(x>7)
	B= stats.linregress(x[mask],y[mask])
	x_in=np.linspace(5,9,100)				
	y_in=B[0]*x_in+B[1]
	print B
	plt.plot(x_in,y_in)
	plt.show()

########################################################################
# verification of the DECRAESE_EMMISIVITY_AFTER_TLIFE flag
########################################################################
def lower(t,E0):
	return E0
def upper(t,tlife,E0):
	return E0*np.power(t/tlife ,-4.)

def verif():
	tlife = 10**(6.565)
	E0 = 10**(52.71)
	
	n = 10000
	t = np.linspace(1e4,1e9,n)
	y = np.zeros(n)

	E0 = 1.
	print integrate.quad(lower, 0, tlife, args=(E0,))
	print integrate.quad(upper, tlife, np.inf,  args=(E0,tlife,))
	
	for i in range(n):
		if t[i]<=tlife :
			y[i]=lower(t[i],E0)
		else:
			y[i]=upper(t[i],tlife,E0)

	plt.plot(np.log10(t),np.log10(y))

	"""
	mask = np.where(t<tlife)
	y[mask] = E0 * np.ones(len(mask))
	
	mask = np.where(t>=tlife)
	y[mask] = E0*np.power(t/tlife,-4)* np.ones(len(mask))
	
	plt.plot(np.log10(t),np.log10(y))
	"""
	
	plt.show()


def verif100000K():
	x=np.linspace(1e2,1e4,1e4)*1e-10
	y=black_body(x,50000)/x *1e7 *1e6

	plt.loglog(x*1e10,y,'r')
	plt.show()
	
	
	
def getNphot(t):
	data = np.loadtxt('python/Starburst99/912_inst_e.dat.txt',unpack=True)
	x=data[0]
	y=np.power(10,data[1])
	return np.interp(t,x,y)	
	
########################################################################
## Atomic.h with starburst99 file
########################################################################


def cross_section(egy):
  P=2.963
  x=egy/4.298e-1
  
  a=(x-1.)**2
  b=np.power(x,P/2.-5.5)
  c=np.power(1+np.sqrt(x/3.288e1),-P)
    
  return 5.475e4*a*b*c*1e-22*(egy >= 13.6)	# m2

def integEnergy():	
	filename='python/Starburst99/fig7e.dat.txt'

	#get constant
	c = Constantes()

	#reading header
	file = open(filename, 'r')
	file.readline()
	file.readline()
	name=file.readline().split()[2:]
	file.close()
	age = [	float(n[:-3]) for n in name]
	
	#loading data
	data = np.loadtxt(filename,unpack=True, skiprows=3)
	
	#converting in SI
	X = data[0]*1e-10 #A->m
	data[1:]=np.power(10,data[1:]) *1e-7/1e-10 #log(Erg/A) -> J/m/s
	
	for num in range(1,len(name)+1):

		#selecting spectrum
		Y=data[num]

		#keep only ionising photons
		mask=np.where(X<=912*1e-10) #912A = 13.6eV
		y=Y[mask]
		x=X[mask]
		
		#integrate total ionising energy
		Eion= integrate.trapz(y,x) #J/s

		#computing mean photon energy
		Nphot_per_sec = getNphot(age[num-1]*1e6)
		#Nphot_by_sec = integrate.trapz(y/lambda2E(x),x)
		Epho= Eion/Nphot_per_sec
		
		#computing nphot/sec/kg
		nphot=Nphot_per_sec/(1e6*c.M0)
	
		#compute photoionisation cross section
		s_lambda = cross_section(lambda2E(x)/c.eV)
		
		#compute sigma N
		N_lambda = y/lambda2E(x)
		s_N = integrate.trapz(s_lambda*N_lambda,x) /Nphot_per_sec

		#compute sigma E
		s_E = integrate.trapz(s_lambda*y,x)/Eion				
		
		#printing
		#print "%3d Myr    Nphot=%e    hnu=%2f*%s    sn=%e    se=%e"%(age[num-1],nphot, Epho/1.602e-19,1.602e-19,s_N,s_E )
		
		print "#ifdef STARBURST_%dMYR"%(age[num-1])
		print "#define SECTION_EFFICACE hnu[0]=%2f*1.6022e-19;alphae[0]=%e*c;alphai[0]=%e*c;"%(Epho/c.eV,s_N,s_E )
		print "#define FACTGRP factgrp[0]=1.0;"
		print "#define SRCINT (%e)"%(nphot)
		print "#endif"
		
		#ploting
#		plt.loglog(X,Y,'k--')
#		plt.loglog(x,y, label="%d Myr"%(age[num-1]))
#		plt.legend()
#	plt.show()
	

