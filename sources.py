
# coding: utf-8

# In[1]:

get_ipython().magic('matplotlib notebook')
import matplotlib.pyplot as plt
import numpy as np
import sys

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import spatial
import time
import scipy.linalg
import healpy as hp

sys.path.insert(0,"/home/deparis/jupyter/pyemma/")
from pyemma import *
get_ipython().magic('cd "~/Emma"')

get_ipython().magic('load_ext autoreload')
get_ipython().magic('autoreload 2')
get_ipython().magic('connect_info')


# In[2]:

runset=db.Runset()
runset.load()
runset.get_description()
runset.get_folder()


# In[327]:

run1=io.Run(runset.runs[4].folder)
run2=io.Run(runset.runs[0].folder)
run3=io.Run(runset.runs[1].folder)
run4=io.Run("/home/deparis/curie_data/data/oct_src/")


# In[109]:

# cur_step =run1.step_00017
# cur_step =run2.step_00021
cur_step =run3.step_00018


# In[110]:

x=cur_step.grid.x.data
y=cur_step.grid.y.data
z=cur_step.grid.z.data
l=cur_step.grid.l.data
d=cur_step.grid.rfield_src.data

projection_level=11
map1=grid.get_cube(x,y,z,l,d,projection_level,"2d")


# In[111]:

plt.figure(figsize=(13,10))
plt.imshow(np.log10(map1),interpolation="nearest",origin="lower")


# In[90]:

# plt.figure()
# bins=np.logspace(-5,0,16)
# x=(bins[1:]+bins[:-1])/2
# dx =np.diff(bins)

# l=run1.step_00017.grid.l.data
# dv = np.power(0.5,3*l)
# y,_=np.histogram(run1.step_00017.grid.rfield_src.data*dv, bins=bins)
# # y=np.cumsum((y*x)[::-1])[::-1]
# y=np.cumsum(x*y)
# plt.plot(x,y,'o:', label="1/8")



# l=run2.step_00021.grid.l.data
# dv = np.power(0.5,3*l)
# y,_=np.histogram(run2.step_00021.grid.rfield_src.data*dv, bins=bins)

# # y=np.cumsum((y*x)[::-1])[::-1]
# y=np.cumsum(x*y)
# plt.plot(x,y,'ro:', label="1/512")

# plt.legend(frameon=False, loc=4 )
# plt.xlabel("src*dv")
# plt.ylabel("#")
# plt.xscale("log")
# # plt.yscale("log")

plt.figure()

mask = run1.step_00017.grid.rfield_src.data !=0
dv = np.power(0.5,3*run1.step_00017.grid.l.data)

d= (run1.step_00017.grid.rfield_src.data*run1.step_00017.grid.field_d.data)[mask]

args=np.argsort(d)
y=np.cumsum(d[args])
x=np.arange(0,1,1./(len(y)))
plt.plot(y,x,'b', label="1/8")


# l=run2.step_00021.grid.l.data
# dv = np.power(0.5,3*l)
# d= run1.step_00017.grid.rfield_src.data*dv
# y=np.cumsum(np.sort(d))
# plt.plot(y, label="1/512")


mask =run2.step_00021.grid.rfield_src.data !=0
dv = np.power(0.5,3*run2.step_00021.grid.l.data)
d= (run2.step_00021.grid.rfield_src.data*run2.step_00021.grid.field_d.data)[mask]

args=np.argsort(d)
y=np.cumsum(d[args])
x=np.arange(0,1,1./(len(y)))
plt.plot(y,x,'r', label="1/512")

plt.legend(frameon=False, loc=4 )
plt.xlabel("cumsum(src*dv) ")
plt.ylabel("n/ntot")
plt.xscale("log")
# plt.ylim(0.5,1.1)


# In[412]:

plt.figure()



# mask =run2.step_00021.grid.rfield_src.data !=0
# dv = np.power(0.5,3*run2.step_00021.grid.l.data[mask])
# x=run2.step_00021.grid.field_d.data[mask]
# y=run2.step_00021.grid.rfield_src.data[mask]*dv
# plt.plot(x,y,'r.', label="1/512")

mask = run1.step_00017.grid.rfield_src.data !=0
dv = np.power(0.5,3*run1.step_00017.grid.l.data[mask])
x=run1.step_00017.grid.field_d.data[mask]
y=run1.step_00017.grid.rfield_src.data[mask]*dv
plt.plot(x,y,'b.', label="1/8")



# mask =run3.step_00018.grid.rfield_src.data !=0
# dv = np.power(0.5,3*run3.step_00018.grid.l.data[mask])
# x=run3.step_00018.grid.field_d.data[mask]
# y=run3.step_00018.grid.rfield_src.data[mask]*dv
# plt.plot(x,y,'g.', label="1/8 oct")



run=run1

#seuil
plt.axvline(run.param.info.ob/run.param.info.om *50, c="k", ls="--")

#emmisivité d'une etoile jeune
mstar=run.param.info.mass_res_star*1.9891e30
src = mstar *run.param.run.src_int_or_fesc
src1star = src/ np.power(run.param.info.unit_l, 3)  *run.param.info.unit_t *0.15**2 /1.605487
plt.axhline(src1star, c="k", ls="--")

#rho^1.5
x=np.logspace(-2,6,50)
y=1e-9*np.power(x,1.5)
plt.plot(x,y,'k--')

#emmisivité a t=explosion SN
age = run.param.run.tlife_SN
y=src1star * np.power(age/ run.param.run.tlife  ,-4.)
plt.axhline(y, c="k", ls="--")

#emmisivité a t=100 tlife
age = run.param.run.tlife *100
y=src1star * np.power(age/ run.param.run.tlife  ,-4.) *run.param.run.ejecta_proportion
plt.axhline(y, c="k", ls="--")





nBins = 128
# l = np.linspace(0.0001, 5, nBins)
l = np.logspace(-6, 0, nBins)

Ms = run.param.info.mass_res_star*1.9891e30/run.param.info.unit_mass
e  = run.param.run.eff_or_tcar
dv = 0.5**(3*11)
dt = 1.192093e-07

rho1 = (l * Ms / e / dt / dv)**(2./3)

pS = np.zeros(nBins)
for i in range(nBins):
    pS[i] = poisson( l[i], np.arange(20) )[1:].sum()

dt*= run.param.info.unit_t *0.15**2 /(365*24*3600)
plt.plot( rho1, fluxUV(1./pS*dt), 'ro--')


plt.legend(frameon=False, loc=4)
plt.xlabel("rho")
plt.ylabel("src")
plt.xscale("log")
plt.yscale("log")
# plt.ylim(0.5,1.1)


# In[ ]:




# In[411]:

def fluxUV(age):
    y=np.zeros(len(age))
    y[age <= run.param.run.tlife] = src1star         
    y[age >  run.param.run.tlife] = src1star*np.power(age[age > run.param.run.tlife]/run.param.run.tlife  ,-8.)
    return y


# In[368]:

import scipy.misc
def poisson( l, k ):
    return l**k * np.exp(-l) / scipy.misc.factorial(k)


# In[295]:

plt.figure()
for i in range(nBins):
    plt.plot( src1star*k, poisson( l[i], k ), '-o' )
    


# In[365]:

p0 = np.zeros(nBins)
pS = np.zeros(nBins)

for i in range(nBins):
    p0[i] = poisson( l[i], k )[0]
    pS[i] = poisson( l[i], k )[1:].sum()
# pS = 1.-p0

plt.figure()
# plt.plot( l, p0, '-o' )
plt.plot( rho1, fluxUV(1./pS*dt), '-o' )


# In[341]:

dt=1e6#*(365*24*3600)

fluxUV(1./pS*dt)


# In[331]:

dt = 3.834559e-04
print(dt/ (run.param.info.unit_t *0.15))


# In[302]:

plt.figure()
# plt.plot( rho1, p0, '-o' )
plt.plot( rho1, pS, '-o' )
# plt.plot( rho2, p0, '-o' )
plt.plot( rho2, pS, '-o' )
plt.xscale('log')
plt.yscale('log')


# In[276]:

run.param.info.mass_res_star*1.9891e30/run.param.info.unit_mass


# In[277]:

run.param.run.eff_or_tcar


# In[262]:

# cur_step =run1.step_00017
cur_step =run2.step_00021
# cur_step =run3.step_00018


mask = cur_step.grid.rfield_src.data !=0
dv = np.power(0.5,3*cur_step.grid.l.data[mask])
y=np.log10(cur_step.grid.field_d.data[mask])
x=np.log10(cur_step.grid.rfield_src.data[mask]*dv)

h,bx,by=np.histogram2d(x,y,bins=128)

plt.figure()
plt.imshow(np.log10(h),interpolation="none", origin='lower', cmap="hot")

# plt.xscale("log")
# plt.yscale("log")

# cbar=plt.colorbar()
# cbar.set_label(r"log$_{10}( Volume ) $")

# plt.clim(-10,-2)
# plt.xlabel(r"log$_{10} (\delta) $")
# plt.ylabel("z")
# plt.xlim(-2,6)
# plt.ylim(5,16)


# In[4]:

cur_step =run2.step_00021

cur_step.fof.get_star(cur_step.star)
cur_step.fof.get_time_oldest_star(cur_step.star)
cur_step.fof.get_time_newest_star(cur_step.star)


# In[14]:

mask=cur_step.fof.time_newest_star!=0
y=cur_step.fof.part_mass_fine[mask]
# x=t-cur_step.fof.time_oldest_star[mask]
x=t-cur_step.fof.time_newest_star[mask]

plt.figure()
plt.plot(x,y,'.')
# plt.xscale('log')
plt.yscale('log')


# In[21]:

t=time.a2t_quad(0.15, run1.param.info.om, run1.param.info.H0)
t=np.min(cur_step.fof.time_newest_star)

mask=cur_step.fof.time_newest_star>0

x=cur_step.fof.time_newest_star[mask]
y=cur_step.fof.part_mass_fine[mask]

print(np.min(x),np.max(x))
print(np.min(y),np.max(y))

nbins=64
xbins=np.logspace(np.log10(np.min(x)),np.log10(np.max(x)), nbins)
ybins=np.logspace(np.log10(np.min(y)),np.log10(np.max(y)), nbins)

# xbins=np.linspace(np.min(x),np.max(x), nbins)
# ybins=np.linspace(np.min(y),np.max(y), nbins)

h,_,_=np.histogram2d(y,x,bins=(ybins,xbins))

plt.figure()
extent = (np.min(x),np.max(x),np.min(y),np.max(y))
plt.imshow(np.log10(h),interpolation="none", origin='lower', extent=extent)
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Mhalo [Mo]')
plt.xlabel('Age')


# In[442]:

plt.figure()

for run in [run0,run1,run2,run3]:
    x=run.param.avg.z
    y=run.param.avg.SFR
    plt.semilogy(x,y, label=run.param.info.mass_res_star)

plt.legend()
plt.xlim(5,15)


# In[330]:

plt.figure()

# for run in [run2,run3,run0,run1]:
for run in [run2,run4]:
    x=run.param.avg.z
    y=run.param.avg.xion.mean
    plt.semilogy(x,y, label=run.param.info.mass_res_star)

plt.legend()
plt.xlim(5,15)


# In[142]:

cur_step= run2.step_00021
cur_step.fof.get_star(cur_step.star)
cur_step.fof.get_luminosity_1600(cur_step)
cur_step.fof.get_luminosity_UV(cur_step)


# In[143]:

nbins=16

# mag1 = cur_step.fof.mag_1600
mag1 = run1.step_00017.fof.mag_UV
mag2 = run2.step_00021.fof.mag_UV

mag_min= np.min(mag1)
mag_max= np.max(mag1)

bins=np.linspace(mag_min,mag_max,nbins+1)
_x=(bins[1:]+bins[:-1])/2
dx=np.diff(bins)

################################################################################

n1,_=np.histogram(mag1,bins=bins)
n2,_=np.histogram(mag2,bins=bins)


box_V=(run1.param.info.box_size_hm1_Mpc/0.67)**3

fig=plt.figure(figsize=(10,9))
fig.patch.set_facecolor('#ffffff')
fig.patch.set_alpha(1)

plt.plot(_x,n1/dx/box_V ,label="1/8",lw=2)
plt.plot(_x,n2/dx/box_V ,label="1/64",lw=2)

# x,y=observations.luminosity_function_fit(6)
# plt.plot(x,y,'k--', label="z=6 fits")

plt.ylim(1e-4,1e-0)
plt.yscale("log", nonposy='mask')
plt.legend(loc=0, frameon=False)

fig.get_axes()[0].invert_xaxis()

plt.title("z=5.7")
plt.xlabel('Mag')
plt.ylabel('$dN/dMag  [h.cMpc^{-1}]$')


# In[426]:

nbins=16

mask = run1.step_00017.grid.rfield_src.data !=0
mag1 = run1.step_00017.grid.field_d.data[mask]

# mag_min= np.min(mag1)
# mag_max= np.max(mag1)
# bins=np.linspace(mag_min,mag_max,nbins+1)
bins=np.logspace(1,6,nbins+1)
_x=(bins[1:]+bins[:-1])/2
dx=np.diff(bins)

################################################################################

n1,_=np.histogram(mag1,bins=bins)

fig=plt.figure(figsize=(10,9))
plt.plot(_x,n1,'o--')
plt.xscale('log')
plt.yscale('log')
plt.xlabel("rho")
plt.ylabel("#")


# In[448]:

#UV luminosity function of 1600 luminosity

cur_step.fof.get_star(cur_step.star)
cur_step.fof.get_luminosity_1600(cur_step)
cur_step.fof.get_luminosity_UV(cur_step)

plt.figure()
x=cur_step.fof.star_flux_1600
y=cur_step.fof.star_flux_UV
plt.plot(x,y,'.')

plt.xscale('log')
plt.yscale('log')

plt.xlabel("1600")
plt.ylabel("UV")


# $e= \frac{1}{2} \rho V^2 +\frac{3}{2} \rho K T$

# In[159]:

plt.figure()


labels=["1/512","1/64"]
colors=["r","b"]
for i,cur_step in enumerate([run1.step_00021, run3.step_00017]):
    
    d=cur_step.grid.field_d.data
    l=cur_step.grid.l.data
    dv=np.power(0.5,3*l)
    src =cur_step.grid.rfield_src.data

    mask= (src!=0) 
    x=src[mask]*dv[mask]
    
    if i==0:
        nbins=64
        xbins=np.logspace(np.log10(np.min(x)),np.log10(np.max(x)), nbins)
        _x=(xbins[1:]+xbins[:-1])/2
    
    h,_=np.histogram(x, bins=xbins)
# np.cumsum(h[::-1])[::-1]
    plt.plot(_x,h, label=labels[i], c=colors[i])

plt.legend()
plt.xscale('log')
plt.yscale('log')


# In[200]:

x=cur_step.star.idx.data
print(len(x))
print(len(np.unique(x)))
print(x)


# In[322]:

print(run1.param.info.mass_res_star)
print(run4.param.info.mass_res_star)


# In[329]:

f, ax = plt.subplots(1,2, sharey=True, figsize=(10,10))

# for i,cur_step in enumerate([run3.step_00017, run1.step_00021]):
for i,cur_step in enumerate([run4.step_00018, run2.step_00018]):
    d=cur_step.grid.field_d.data
    p=cur_step.grid.field_p.data
    u=cur_step.grid.field_u.data
    v=cur_step.grid.field_v.data
    l=cur_step.grid.l.data
    w=cur_step.grid.field_w.data
    T=cur_step.grid.rfield_temp.data
    src =cur_step.grid.rfield_src.data
    xion =cur_step.grid.xion.data
    dv=np.power(0.5,3*l)

#     ec = 1./2 * d* (np.power(u,2)+np.power(v,2)+np.power(w,2))
#     ei = 3./2 * d* T
#     et = ec+ei

    mask= (src!=0) & (xion!=1)
    # x=src[mask]
    # y=et[mask]

    # x=d[mask]
    # y=src[mask]/et[mask]

    # x=d
    # y=et

    # x=T[mask]*d[mask]*dv[mask]
    # y=src[mask]*dv[mask]
    # # z=dv[mask]
    # z=None

    # x=np.sqrt(T[mask]/d[mask])

    

    
#     x=xion[mask]
    
    
    
    y=src[mask]*dv[mask]
    
    x = d[mask] * np.power(xion[mask],2)/(1.-xion[mask])
    x=xion[mask]/(1-xion[mask])
#     x=T[mask]
#     x=CompCooling(T[mask],xion[mask],d[mask]*(1-xion[mask]))
#     x=p[mask]
    print(x)
    
    nh2=d[mask]*xion[mask]
    temp=T[mask]
    alpha=  1.778e-29 * temp*(2.*157807./temp)**1.965 / (1.+(2.*157807./temp/0.541)**0.502)**2.697 * xion[mask]**2 * nh2**2 

    x=alpha*np.power(nh2,2) * np.power(xion[mask],2)
            
    x=d[mask]
    
    
    z=dv[mask]
    
    x=x[x!=0]
    y=y[x!=0]
    z=z[x!=0]
    

    nbins=100
    xbins=np.logspace(np.log10(np.min(x)),np.log10(np.max(x)), nbins)
    ybins=np.logspace(np.log10(np.min(y)),np.log10(np.max(y)), nbins)
    extent = (np.min(x),np.max(x),np.min(y),np.max(y))

    h,bx,by=np.histogram2d(y,x,bins=(ybins,xbins),weights=z)
    ax[i].imshow(np.log10(h),interpolation="none", origin='lower', cmap="hot", extent=extent, aspect='auto')

    ax[i].set_xscale('log')
    ax[i].set_yscale('log')
    
    
    if i==0:
        xmin=np.min(x)
        xmax=np.max(x)
        ymin=np.min(y)
        ymax=np.max(y)
    else:
        xmin=min(xmin,np.min(x))
        xmax=max(xmax,np.max(x))
        ymin=min(ymin,np.min(y))
        ymax=max(ymax,np.max(y))

for i in range(len(ax)):
    ax[i].set_xlim(xmin,xmax)
    ax[i].set_ylim(ymin,ymax)

ax[0].set_ylabel("src.dv")
for i in range(len(ax)):
    ax[i].set_xlabel("x/(1-x)")
    ax[i].set_xlabel("recombination rate")
#     ax[i].set_xlabel("T")
                 
f.subplots_adjust(wspace=0)
# plt.setp([a.get_yticklabels() for a in f.axes[:-1]], visible=False)


# In[294]:

def CompCooling( temp, x, nH, aexp=1., CLUMPF=1. ):
    """
    return:
       lambda_cool [J.m^-3.s^-1]
       tcool [Myr]
    """

    
    nh2 = nH*1e-6 ### unit convertion [m^-3] ==> [cm^-3]
    
    ### Collisional Ionization Cooling
    c1 = 1.27e-21 * np.exp(-157809./temp) * np.sqrt(temp)/(1.+np.sqrt(temp/1e5)) * x*(1.-x) * nh2**2 *CLUMPF
    
    ### Case A Recombination Cooling
    c2 = 1.778e-29 * temp*(2.*157807./temp)**1.965 / (1.+(2.*157807./temp/0.541)**0.502)**2.697 * x**2 * nh2**2 *CLUMPF
    
    ### Case B Recombination Cooling
    c6 = 3.435e-30 * temp*(2.*157807./temp)**1.970 / (1.+(2.*157807./temp/2.250)**0.376)**3.720 * x**2 * nh2**2 *CLUMPF
    
    ### Collisional excitation cooling
    c3 = np.exp(-118348./temp)*7.5e-19/(1+np.sqrt(temp/1e5))*x*(1.-x)*nh2**2 *CLUMPF
    
    ### Bremmsstrahlung
    c4 = 1.42e-27 * 1.5 * np.sqrt(temp) * x**2 * nh2**2 *CLUMPF
    
    ### Compton Cooling
    c5 = 1.017e-37 * (2.727/aexp)**4 * (temp-2.727/aexp)* nh2*x
    
    ### Overall Cooling
    lambda_cool = c1+c2+c3+c4+c5+c6 ### [erg.cm^-3.s^-1]
    lambda_cool = lambda_cool * 1e-7 *1e6 ### unit convertion [erg.cm^-3.s^-1] ==> [J.m^-3.s^-1]
    ### cooling times
    unsurtc = np.amax( (c1, c2, c3, c4, np.abs(c5), c6), axis=1 ) * 1e-7 ### [J.cm^-3.s^-1]
    unsurtc=c6
    tcool = 1.5 * nh2 * (1+x) * 1.38064852e-23 * temp/unsurtc ### [Myr]
    return lambda_cool#, tcool


# In[283]:

get_ipython().magic('pinfo np.max')

