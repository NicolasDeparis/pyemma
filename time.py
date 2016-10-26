# coding: utf-8

import numpy as np
from scipy import integrate

def a2t_quad(az, o_m, H0):

    o_v=1.-o_m

    H0 = H0*1e3/1e6/3.08567758e16*(365*24*3600) # Hubble constant in SI unit

    def f(a):
        return a/ np.sqrt(o_m*a + o_v*a**4 )

    return 1./H0 * integrate.quad(f,0,az)[0]

def a2t_romberg(b,info):


	MAXJ=5
	MAXITER=16

	H0=info.H0 *(365*24*3600)
	om=info.om
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

	g=np.empty(MAXJ+2)
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
