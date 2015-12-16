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
		self.yr = 31556926						# secondes in 1 year
		
		
		
		
		#cosmo
		
		self.H0 = 67.0			     			# Hubble constant 
		self.H0_SI = self.H0*1e3/1e6/self.Parsec # Hubble constant in SI unit
		self.h  = self.H0/100.
	
		self.WB = 0.049
		self.WM = 0.3175                 		# Omega(matter)
		self.WV = 0.6825               			# Omega(vacuum) or lambda
		self.WR = 0#1*4.165E-5/(self.h*self.h)    	# Omega(radiation)
		self.WK = 0#1-self.WM-self.WR-self.WV		# Omega curvaturve = 1-Omega(total)

		


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


def a2t(az):
	""" 
#		convert expansion factor to time
	"""
	c = Constantes()
	
	age = 0.
	n=1000         # number of points in integrals
	for i in range(n):
		a = az*(i+0.5)/n
		adot = np.sqrt(c.WK+(c.WM / a)+(c.WR/ (a*a) )+ (c.WV*a*a) )
		age = age + 1./adot
	zage = az*age/n
	zage = (977.8 /c.H0)*zage * 1e9 

	return zage


def a2t_romberg(b):
	
	c = Constantes()
	
	MAXJ=5
	MAXITER=16
	
	H0=c.H0_SI *(365*24*3600)
	om=c.WM	
	ov=1.-om	
	a=1e-8
	tol=1e-8
  	
	def faexp(a, om, ov):
		return a*a/np.sqrt((1.0-om-ov)*a*a*a*a+om*a*a*a+ov*a*a*a*a*a*a)	
	
	def JMAX(i,MAXJ):		
		if i<MAXJ:
			return i
		else:
			return MAXJ		

	g=np.zeros(MAXJ+2)
	h=0.5*(b-a)	
	gmax=h*(faexp(a,om,ov)+faexp(b,om,ov))
	g[1]=gmax;
	nint=1;
	error=np.inf;
	

	i=0;
	while True:
		i+=1
		if( (i>MAXITER) or ((i>5) and (np.abs(error)<tol))): 
			break
		
		g0=0.
		for k in range(1,nint+1):		
			g0+=faexp(a+(k+k-1)*h,om,ov)

		g0=0.5*g[1]+h*g0
		h*=0.5
		nint=nint+nint
		
		jmax=JMAX(i,MAXJ)
				
		fourj=1.0
		for j in range(1,jmax+1):		
			fourj*=4.0
			g1=g0+(g0-g[j])/(fourj-1.0)
			g[j]=g0
			g0=g1
		
		if(np.abs(g0)>tol):
			error=1.0-gmax/g0
		else:
			error=gmax
		
		gmax=g0
		g[jmax+1]=g0

	res=g0

	if((i>MAXITER) and (np.abs(error)>tol)):
		print("ROMBINT FAILED TO CONVERGE integ=%e error=%e\n"%(res,error) )
  
	return res/H0




def a2t_quad(az):
	c = Constantes()
	
	o_m=c.WM	
	o_v=1.-o_m
	H0=c.H0_SI *(365*24*3600)
	
	def f(a):
		return a/ np.sqrt(o_m*a + o_v*a**4 )
		
	return 1./H0 * integrate.quad(f,0,az)[0]


def t2a(t):
	"""
		convert time to expansion factor
	"""
	n = 10000
	A = np.arange(n+1) / float(n)
	T = a2t(A)
	return np.interp(t,T,A)
	
def code2t(step,t):
	"""
		convert code time into physical time
	"""
	p =step.param.info.get()
	a= step.a
	scale = 1./(a**2 *p["unit_t"])
	
	return t/scale
		

########################################################################
# density
########################################################################

def code2d(step,rho):
	"""
		convert code density to physical density
	"""
	p =step.param.info.get()
		
	a= step.a
	
	unit_d = p["unit_mass"]/pow(p["unit_l"],3)
	
	return rho/pow(a,3)*unit_d
		
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
	

def miralda_escude_fit(d, redshift):
    if redshift == 2:
        a=0.406
        c=0.558
        d0=2.54
        b=2.23
    elif redshift == 3:
        a=0.558
        c=0.599
        d0=1.89
        b=2.35
    elif redshift == 4:
        a=0.711
        c=0.611
        d0=1.53
        b=2.48
    elif redshift == 6:
        a=0.864
        c=0.880
        d0=1.09
        b=2.50
    elif redshift == -6:
        #Pawlik fit at z=6 from http://arxiv.org/pdf/0807.3963.pdf
        a = 3.038
        c = -0.932
        d0 = 1.477 
        b = 3.380 
    else:
        print "redshift invalide, available are 2,3,4,6"
        return None
    
    d23 = np.power(d,-2./3.)
    up = np.power(d23-c,2.)
    down = 2. * np.power(2./3.*d0,2.)    
    return a * np.exp(-up/down) *np.power(d,-b)
