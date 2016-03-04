
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


# In[4]:

runset=db.Runset()
runset.load()
runset.get_description()
runset.get_folder()


# In[23]:

run1=io.Run(runset.runs[10].folder) 


# In[24]:

print(runset.runs[0].labels)
cur_run=run1
cur_step=cur_run.step_00017


# In[17]:

run2=io.Run(runset.runs[4].folder)


# In[18]:

print(runset.runs[1].labels)
cur_run=run2
cur_step=cur_run.step_00016


# In[27]:

run3=io.Run(runset.runs[4].folder)


# In[28]:


cur_run=run3
cur_step=cur_run.step_00020


# In[10]:

#cur_cat.gen(8,cur_run.param.info.nproc)

cur_step.fof.read_masst()
cur_step.fof.read_struct()

cur_step.fof.get_R200(cur_run.param.info.ob, cur_run.param.info.om)

#________________________________________________________________________
cur_step.fof.get_part(cur_step.part, force=0)
cur_step.fof.get_part_mass(cur_step.part,cur_run.param.info)

# cur_step.fof.get_part_fine(cur_step.part, force=0)
cur_step.fof.get_part_mass_fine(cur_run.param.info)

#________________________________________________________________________
cur_step.fof.get_star(cur_step.star, force=0)
cur_step.fof.get_star_mass(cur_step.star,cur_run.param.info)

cur_step.fof.get_stars_fine(cur_step.grid, cur_step.star)
cur_step.fof.get_star_mass_fine(cur_step.star,cur_run.param.info)

#________________________________________________________________________
cur_step.fof.get_cells(cur_step.grid, force=0)
cur_step.fof.get_gas_mass(cur_step.grid,cur_run.param.info)

cur_step.fof.get_cells_fine(cur_step.grid)
cur_step.fof.get_gas_mass_fine(cur_step.grid,cur_run.param.info)

#________________________________________________________________________

cur_step.fof.get_mean_vel()


# # DM PROFILE

# In[ ]:

n=3

nbins=32
plt.figure()

################################################################################

arg= np.argsort(cur_step.fof.npart)
halo_num=arg[-n]

xc=cur_step.fof.x[halo_num]
yc=cur_step.fof.y[halo_num]
zc=cur_step.fof.z[halo_num]

cur_step.part.x.read()
cur_step.part.y.read()
cur_step.part.z.read()
cur_step.part.mass.read()

part=cur_step.fof.part[halo_num]
x=cur_step.part.x.data[part]
y=cur_step.part.y.data[part]
z=cur_step.part.z.data[part]
m=cur_step.part.mass.data[part]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

bins=np.linspace(0,np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

V = 4./3.*np.pi * np.power(bins,3)
dV = np.diff(V)

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

# plt.errorbar(_x,h11/dV,yerr=err1)
plt.errorbar(_x,h11/dV,yerr=err1)

plt.yscale("log")
plt.xlabel("Radius [box unit]")
plt.ylabel("DM density")


# # STELLAR PROFILE

# In[8]:

n=2
nbins=32
plt.figure()

################################################################################

arg= np.argsort(cur_step.fof.npart)
halo_num=arg[-n]

cur_step.fof.get_star(cur_step.star)

xc=cur_step.fof.x[halo_num]
yc=cur_step.fof.y[halo_num]
zc=cur_step.fof.z[halo_num]

cur_step.star.x.read()
cur_step.star.y.read()
cur_step.star.z.read()
cur_step.star.mass.read()

stars=cur_step.fof.stars[halo_num]
x=cur_step.star.x.data[stars]
y=cur_step.star.y.data[stars]
z=cur_step.star.z.data[stars]
m=cur_step.star.mass.data[stars]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

V = 4./3.*np.pi * np.power(bins,3)
dV = np.diff(V)

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

plt.errorbar(_x,h11/dV,yerr=err1)

plt.legend()
plt.yscale("log")
plt.xlabel("Radius [box unit]")
plt.ylabel("Stellar density")


# In[ ]:

n=4
nbins=32
plt.figure()

################################################################################

arg= np.argsort(cur_step.fof.npart)
halo_num=arg[-n]

xc=cur_step.fof.x[halo_num]
yc=cur_step.fof.y[halo_num]
zc=cur_step.fof.z[halo_num]

cur_step.grid.x.read()
cur_step.grid.y.read()
cur_step.grid.z.read()
cur_step.grid.l.read()
cur_step.grid.field_d.read()

cells=cur_step.fof.cells[halo_num]
x=cur_step.grid.x.data[cells]-xc
y=cur_step.grid.y.data[cells]-yc
z=cur_step.grid.z.data[cells]-zc
l=cur_step.grid.l.data[cells]
dv=np.power(0.5,3*l)
m=cur_step.grid.field_d.data[cells]*dv

r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2))

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

V = 4./3.*np.pi * np.power(bins,3)
dV = np.diff(V)

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

# plt.errorbar(_x,h11,yerr=err1)
plt.errorbar(_x,h11/dV,yerr=err1)

plt.legend()
plt.yscale("log")
plt.xlabel("Radius [box unit]")
plt.ylabel("Gas density")


# # HYDRO FLUX PROFILE

# In[ ]:

n=4
nbins=32
plt.figure()

################################################################################

arg= np.argsort(cur_step.fof.npart)
halo_num=arg[-n]

xc=cur_step.fof.x[halo_num]
yc=cur_step.fof.y[halo_num]
zc=cur_step.fof.z[halo_num]

cur_step.grid.x.read()
cur_step.grid.y.read()
cur_step.grid.z.read()
cur_step.grid.field_d.read()
cur_step.grid.field_u.read()
cur_step.grid.field_v.read()
cur_step.grid.field_w.read()

cells=cur_step.fof.cells[halo_num]
x=cur_step.grid.x.data[cells]-xc
y=cur_step.grid.y.data[cells]-yc
z=cur_step.grid.z.data[cells]-zc
d=cur_step.grid.field_d.data[cells]
u=cur_step.grid.field_u.data[cells]
v=cur_step.grid.field_v.data[cells]
w=cur_step.grid.field_w.data[cells]

r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2))

m=(x*u+y*v+z*w)/r 
m*= 4*np.pi*np.power(r,2)

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

plt.errorbar(_x,-h11,yerr=err1)

plt.legend()
#plt.yscale("log")
plt.xlabel("Radius [box unit]")
plt.ylabel("Radial Inflow")
plt.title("Hydrodynamical radial Inflow")


# # RADIATION OUTFLOW

# In[19]:

n=2


nbins=32

################################################################################

arg= np.argsort(cur_step.fof.npart)
halo_num=arg[-n]

xc=cur_step.fof.x[halo_num]
yc=cur_step.fof.y[halo_num]
zc=cur_step.fof.z[halo_num]

cur_step.grid.x.read()
cur_step.grid.y.read()
cur_step.grid.z.read()
cur_step.grid.rfield_fx0.read()
cur_step.grid.rfield_fy0.read()
cur_step.grid.rfield_fz0.read()

cells=cur_step.fof.cells[halo_num]
x=cur_step.grid.x.data[cells]-xc
y=cur_step.grid.y.data[cells]-yc
z=cur_step.grid.z.data[cells]-zc
u=cur_step.grid.rfield_fx0.data[cells]
v=cur_step.grid.rfield_fy0.data[cells]
w=cur_step.grid.rfield_fz0.data[cells]

r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2))

m=(x*u+y*v+z*w)/r 
m*= 4*np.pi*np.power(r,2)

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2


n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

plt.figure()
plt.errorbar(_x,h11,yerr=err1)

#plt.yscale("log")
plt.ylabel("Radial outflow")
plt.title("Radiation radial outflow")


# f, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(10,10))

# ax1.errorbar(_x, h11)
# # ax1.set_yscale('log')
# ax1.set_ylabel("Outflow")

# ax2.errorbar(_x, -h11)
# # ax2.set_yscale('log')
# ax2.set_ylabel("Inflow")

# # max_bound= max(np.max(h11[h11>0]),-np.min(h11[h11<0]))
# # min_bound= 10**np.int(np.log10(min(np.min(h11[h11>0]),-np.max(h11[h11<0]))))

# # ax1.set_ylim(min_bound, max_bound)
# # ax2.set_ylim(min_bound, max_bound)
# ax2.invert_yaxis()



# f.subplots_adjust(hspace=0)
# plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
# plt.xlabel("Radius [box unit]")


# # HYDRO

# In[154]:

cur_step.grid.x.read()
cur_step.grid.y.read()
cur_step.grid.z.read()
cur_step.grid.field_d.read()
cur_step.grid.field_u.read()
cur_step.grid.field_v.read()
cur_step.grid.field_w.read()

Nbins=3

mh=cur_step.fof.part_mass_fine
Mbins=np.logspace(np.log10(np.min(mh)), np.log10(np.max(mh)), Nbins+1)

r_tot=np.empty(Nbins,dtype=np.object)
m_tot=np.empty(Nbins,dtype=np.object)

for i in range(Nbins):

    r_bin=[]
    m_bin=[]
    
    for halo_num in np.where( (mh>=Mbins[i])  & (mh<Mbins[i+1]) )[0] :
        
        xc=cur_step.fof.x[halo_num]
        yc=cur_step.fof.y[halo_num]
        zc=cur_step.fof.z[halo_num]

        cells=cur_step.fof.cells[halo_num]
        x=cur_step.grid.x.data[cells]-xc
        y=cur_step.grid.y.data[cells]-yc
        z=cur_step.grid.z.data[cells]-zc
        
        d=cur_step.grid.field_d.data[cells]
        u=cur_step.grid.field_u.data[cells]*d
        v=cur_step.grid.field_v.data[cells]*d
        w=cur_step.grid.field_w.data[cells]*d

        r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2))  / cur_step.fof.R200[halo_num]
        m=(x*u+y*v+z*w)/r * (4*np.pi*np.power(r,2))
                
        r_bin.append(r)
        m_bin.append(m)

    r_tot[i]=np.concatenate(r_bin)
    m_tot[i]=np.concatenate(m_bin)


# In[155]:

plt.figure()

for i in range(Nbins):

    nbins=16
    bins=np.linspace(np.min(r_tot[i]),np.max(r_tot[i]),nbins+1)
    _x=(bins[1:]+bins[:-1])/2

    n1,_=np.histogram(r_tot[i],bins=bins)
    h11,_=np.histogram(r_tot[i],bins=bins,weights=m_tot[i])
    h12,_=np.histogram(r_tot[i],bins=bins,weights=m_tot[i]*m_tot[i])
    err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

    lab= ("%.1e <Mh< %.1e"%(Mbins[i],Mbins[i+1] ))
    plt.errorbar(_x,-h11,yerr=err1, label=lab)

plt.yscale("log")
plt.ylabel("Radial outflow")
plt.title("Hydro radial inflow")
plt.legend(loc=4)
plt.xlim(0,1)


# # RAD

# In[156]:

cur_step.grid.x.read()
cur_step.grid.y.read()
cur_step.grid.z.read()
cur_step.grid.rfield_fx0.read()
cur_step.grid.rfield_fy0.read()
cur_step.grid.rfield_fz0.read()

Nbins=3

mh=cur_step.fof.part_mass_fine
Mbins=np.logspace(np.log10(np.min(mh)), np.log10(np.max(mh)), Nbins+1)

r_tot=np.empty(Nbins,dtype=np.object)
m_tot=np.empty(Nbins,dtype=np.object)

for i in range(Nbins):

    r_bin=[]
    m_bin=[]
    
    for halo_num in np.where( (mh>=Mbins[i])  & (mh<Mbins[i+1]) )[0] :
        
        xc=cur_step.fof.x[halo_num]
        yc=cur_step.fof.y[halo_num]
        zc=cur_step.fof.z[halo_num]

        cells=cur_step.fof.cells[halo_num]
        x=cur_step.grid.x.data[cells]-xc
        y=cur_step.grid.y.data[cells]-yc
        z=cur_step.grid.z.data[cells]-zc
        u=cur_step.grid.rfield_fx0.data[cells]
        v=cur_step.grid.rfield_fy0.data[cells]
        w=cur_step.grid.rfield_fz0.data[cells]

        r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2)) / cur_step.fof.R200[halo_num]
        m=(x*u+y*v+z*w)/r * (4*np.pi*np.power(r,2))
                
        r_bin.append(r)
        m_bin.append(m)

    r_tot[i]=np.concatenate(r_bin)
    m_tot[i]=np.concatenate(m_bin)


# In[157]:

plt.figure()

for i in range(Nbins):

    nbins=16
    bins=np.linspace(np.min(r_tot[i]),np.max(r_tot[i]),nbins+1)
    _x=(bins[1:]+bins[:-1])/2

    n1,_=np.histogram(r_tot[i],bins=bins)
    h11,_=np.histogram(r_tot[i],bins=bins,weights=m_tot[i])
    h12,_=np.histogram(r_tot[i],bins=bins,weights=m_tot[i]*m_tot[i])
    err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

    lab= ("%.1e <Mh< %.1e"%(Mbins[i],Mbins[i+1] ))
    plt.errorbar(_x,h11,yerr=err1, label=lab)

plt.yscale("log")
plt.ylabel("Radial outflow")
plt.title("Radiation radial outflow")
plt.legend(loc=4)

plt.xlim(0,1)


# In[225]:

cur_step.grid.x.read()
cur_step.grid.y.read()
cur_step.grid.z.read()
cur_step.star.x.read()
cur_step.star.y.read()
cur_step.star.z.read()

cur_step.fof.get_part_mass_fine()
cur_step.fof.get_star(cur_step.star)
cur_step.fof.get_luminosity(cur_step)

Nbins=4

mh=cur_step.fof.part_mass_fine
Mbins=np.logspace(np.log10(np.min(mh)), np.log10(np.max(mh)), Nbins+1)
Mbins=np.logspace(7,7+Nbins,Nbins+1)

r_tot=np.empty(Nbins,dtype=np.object)
m_tot=np.empty(Nbins,dtype=np.object)

for i in range(Nbins):

    r_bin=[]
    m_bin=[]
    
    for halo_num in np.where( (mh>=Mbins[i])  & (mh<Mbins[i+1]) )[0] :
        
        xc=cur_step.fof.x[halo_num]
        yc=cur_step.fof.y[halo_num]
        zc=cur_step.fof.z[halo_num]
        stars = cur_step.fof.stars[halo_num]
        
        x=cur_step.star.x.data[stars]-xc
        y=cur_step.star.y.data[stars]-yc
        z=cur_step.star.z.data[stars]-zc
        
        r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2)) / cur_step.fof.R200[halo_num]
        m= cur_step.star.flux[stars] 
                
        r_bin.append(r)
        m_bin.append(m)

    r_tot[i]=np.concatenate(r_bin)
    m_tot[i]=np.concatenate(m_bin)

#     args = np.argsort(r_tot[i])
#     r_tot[i] = (r_tot[i])[args]
#     m_tot[i] = np.cumsum((m_tot[i])[args])


# In[226]:

plt.figure()

for i in range(Nbins):

    nbins=16
#     bins=np.linspace(np.min(r_tot[i]),np.max(r_tot[i]),nbins+1)
    bins=np.linspace(0,1,nbins+1)
    _x=(bins[1:]+bins[:-1])/2

    
    n1,_=np.histogram(r_tot[i],bins=bins)
    h11,_=np.histogram(r_tot[i],bins=bins,weights=m_tot[i])
    h12,_=np.histogram(r_tot[i],bins=bins,weights=m_tot[i]*m_tot[i])
    err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

    h11=np.cumsum(h11)
    lab= ("%.1e <Mh< %.1e"%(Mbins[i],Mbins[i+1] ))
    plt.errorbar(_x,h11,yerr=err1, label=lab)

plt.yscale("log")
plt.ylabel("Stellar flux")
plt.title("")
plt.legend(loc=4)

plt.xlim(0,1)


# In[29]:

cur_step.grid.x.read()
cur_step.grid.y.read()
cur_step.grid.z.read()
cur_step.grid.rfield_fx0.read()
cur_step.grid.rfield_fy0.read()
cur_step.grid.rfield_fz0.read()
cur_step.star.x.read()
cur_step.star.y.read()
cur_step.star.z.read

cur_step.fof.get_cells(cur_step.grid)
cur_step.fof.get_part_mass_fine()
cur_step.fof.get_star(cur_step.star)
cur_step.fof.get_luminosity_UV(cur_step)


Nbins=4

mh=cur_step.fof.part_mass_fine
# Mbins=np.logspace(np.log10(np.min(mh)), np.log10(np.max(mh)), Nbins+1)
Mbins=np.logspace(8,8+Nbins,Nbins+1)

flux_r_tot=np.empty(Nbins,dtype=np.object)
flux_m_tot=np.empty(Nbins,dtype=np.object)

star_r_tot=np.empty(Nbins,dtype=np.object)
star_m_tot=np.empty(Nbins,dtype=np.object)

for i in range(Nbins):

    flux_r_bin=[]
    flux_m_bin=[]
    star_r_bin=[]
    star_m_bin=[]
    
    for halo_num in np.where( (mh>=Mbins[i])  & (mh<Mbins[i+1]) )[0] :
        
        xc=cur_step.fof.x[halo_num]
        yc=cur_step.fof.y[halo_num]
        zc=cur_step.fof.z[halo_num]

        cells=cur_step.fof.cells[halo_num]
        
        x=cur_step.grid.x.data[cells]-xc
        y=cur_step.grid.y.data[cells]-yc
        z=cur_step.grid.z.data[cells]-zc
        
        u=cur_step.grid.rfield_fx0.data[cells]
        v=cur_step.grid.rfield_fy0.data[cells]
        w=cur_step.grid.rfield_fz0.data[cells]

        flux_r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2)) / cur_step.fof.R200[halo_num]
        flux_m=(x*u+y*v+z*w)/flux_r * (4*np.pi*np.power(flux_r,2))
                
        flux_r_bin.append(flux_r)
        flux_m_bin.append(flux_m)
        
        stars=cur_step.fof.stars[halo_num]
        
        x=cur_step.star.x.data[stars]-xc
        y=cur_step.star.y.data[stars]-yc
        z=cur_step.star.z.data[stars]-zc
        
        star_r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2)) / cur_step.fof.R200[halo_num]
        star_m= cur_step.star.flux_UV[stars] 
        
        star_r_bin.append(star_r)
        star_m_bin.append(star_m)
        
    star_r_tot[i]=np.concatenate(star_r_bin)
    star_m_tot[i]=np.concatenate(star_m_bin)
    
    flux_r_tot[i]=np.concatenate(flux_r_bin)
    flux_m_tot[i]=np.concatenate(flux_m_bin)


# In[116]:

plt.figure()

labels=["1/64", "1/8","1/512"]


for istep, cur_step in enumerate([run2.step_00016, run1.step_00017, run3.step_00020]):

#     cur_step.grid.x.read()
#     cur_step.grid.y.read()
#     cur_step.grid.z.read()
#     cur_step.grid.rfield_fx0.read()
#     cur_step.grid.rfield_fy0.read()
#     cur_step.grid.rfield_fz0.read()
#     cur_step.star.x.read()
#     cur_step.star.y.read()
#     cur_step.star.z.read

#     cur_step.fof.get_cells(cur_step.grid)
#     cur_step.fof.get_part_mass_fine()
#     cur_step.fof.get_star(cur_step.star)
#     cur_step.fof.get_luminosity_UV(cur_step)

    Nbins=2

    mh=cur_step.fof.part_mass_fine
    # Mbins=np.logspace(np.log10(np.min(mh)), np.log10(np.max(mh)), Nbins+1)
    Mbins=np.logspace(8,12,Nbins+1)
    
    flux_r_tot=np.empty(Nbins,dtype=np.object)
    flux_m_tot=np.empty(Nbins,dtype=np.object)

    star_r_tot=np.empty(Nbins,dtype=np.object)
    star_m_tot=np.empty(Nbins,dtype=np.object)

    for i in range(Nbins):

        flux_r_bin=[]
        flux_m_bin=[]
        star_r_bin=[]
        star_m_bin=[]

        for halo_num in np.where( (mh>=Mbins[i])  & (mh<Mbins[i+1]) )[0] :

            xc=cur_step.fof.x[halo_num]
            yc=cur_step.fof.y[halo_num]
            zc=cur_step.fof.z[halo_num]

            cells=cur_step.fof.cells[halo_num]

            x=cur_step.grid.x.data[cells]-xc
            y=cur_step.grid.y.data[cells]-yc
            z=cur_step.grid.z.data[cells]-zc

            u=cur_step.grid.rfield_fx0.data[cells]
            v=cur_step.grid.rfield_fy0.data[cells]
            w=cur_step.grid.rfield_fz0.data[cells]

            flux_r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2)) / cur_step.fof.R200[halo_num]
            flux_m=(x*u+y*v+z*w)/flux_r * (4*np.pi*np.power(flux_r,2))

            flux_r_bin.append(flux_r)
            flux_m_bin.append(flux_m)

            stars=cur_step.fof.stars[halo_num]

            x=cur_step.star.x.data[stars]-xc
            y=cur_step.star.y.data[stars]-yc
            z=cur_step.star.z.data[stars]-zc

            star_r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2)) / cur_step.fof.R200[halo_num]
            star_m= cur_step.star.flux_UV[stars] 

            star_r_bin.append(star_r)
            star_m_bin.append(star_m)

        star_r_tot[i]=np.concatenate(star_r_bin)
        star_m_tot[i]=np.concatenate(star_m_bin)

        flux_r_tot[i]=np.concatenate(flux_r_bin)
        flux_m_tot[i]=np.concatenate(flux_m_bin)
    
    color = [ 'b', 'g',  'r']
    for i in range(Nbins):

        nbins=2
        bins=np.linspace(0,1,nbins+1)
        _x=(bins[1:]+bins[:-1])/2

        star_n1,_=np.histogram(star_r_tot[i],bins=bins)
        star_h11,_=np.histogram(star_r_tot[i],bins=bins,weights=star_m_tot[i])
        star_h11/=star_n1

        flux_n1,_=np.histogram(flux_r_tot[i],bins=bins)
        flux_h11,_=np.histogram(flux_r_tot[i],bins=bins,weights=flux_m_tot[i])
        flux_h11/=flux_n1

        star_h11=np.cumsum(star_h11)

        #lab= ("%.1e <Mh< %.1e"%(Mbins[i],Mbins[i+1] ))
        plt.plot(_x,flux_h11/star_h11, ':o',  label=lab, c=color[istep])

                
plt.yscale("log")
plt.ylabel("")


plt.xlim(0,1)
# plt.ylim(1e-57,1e-51)


# In[168]:

plt.figure()

labels=["1/64", "1/8","1/512"]


for istep, cur_step in enumerate([run2.step_00016, run1.step_00017, run3.step_00020]):
# for istep, cur_step in enumerate([run2.step_00016]):

#     cur_step.grid.x.read()
#     cur_step.grid.y.read()
#     cur_step.grid.z.read()
#     cur_step.grid.rfield_fx0.read()
#     cur_step.grid.rfield_fy0.read()
#     cur_step.grid.rfield_fz0.read()
#     cur_step.star.x.read()
#     cur_step.star.y.read()
#     cur_step.star.z.read

#     cur_step.fof.get_cells(cur_step.grid)
#     cur_step.fof.get_part_mass_fine()
#     cur_step.fof.get_star(cur_step.star)
#     cur_step.fof.get_luminosity_UV(cur_step)

    Nbins=2

    mh=cur_step.fof.part_mass_fine
#     Mbins=np.logspace(np.log10(np.min(mh)), np.log10(np.max(mh)), Nbins+1)
    Mbins=np.logspace(8,12,Nbins+1)

    for i in range(Nbins):

        flux_r_bin=[]
        flux_m_bin=[]
        star_r_bin=[]
        star_m_bin=[]
        fesc_r_bin=[]
        fesc_m_bin=[]
        for halo_num in np.where( (mh>=Mbins[i])  & (mh<Mbins[i+1]) )[0] :

            xc=cur_step.fof.x[halo_num]
            yc=cur_step.fof.y[halo_num]
            zc=cur_step.fof.z[halo_num]
#_________________________________________________________________________________________________
            
            cells=cur_step.fof.cells[halo_num]

            x=cur_step.grid.x.data[cells]-xc
            y=cur_step.grid.y.data[cells]-yc
            z=cur_step.grid.z.data[cells]-zc                        
            
            u=cur_step.grid.rfield_fx0.data[cells]
            v=cur_step.grid.rfield_fy0.data[cells]
            w=cur_step.grid.rfield_fz0.data[cells]
            
            r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2)) / cur_step.fof.R200[halo_num]
        
            flux_m=(x*u+y*v+z*w)/r 
            flux_m *= (4*np.pi*np.power(r,2))
            
            nbins=5
            bins=np.linspace(0,1,nbins+1)
            _x=(bins[1:]+bins[:-1])/2

            flux_n1,_=np.histogram(r,bins=bins)
            flux_h11,_=np.histogram(r,bins=bins,weights=flux_m)
            
    #             flux_h11 *= 4**(istep)
            
#_________________________________________________________________________________________________

            stars=cur_step.fof.stars[halo_num]

            x=cur_step.star.x.data[stars]-xc
            y=cur_step.star.y.data[stars]-yc
            z=cur_step.star.z.data[stars]-zc
            r= np.sqrt( np.power(x,2)+np.power(y,2)+np.power(z,2)) / cur_step.fof.R200[halo_num]

            star_m= cur_step.star.flux_UV[stars]          

            star_n1,_=np.histogram(r,bins=bins)
            star_h11,_=np.histogram(r,bins=bins,weights=star_m)

            star_h11=np.cumsum(star_h11)

#_________________________________________________________________________________________________
            
            flux_m_bin.append(flux_h11)
            star_m_bin.append(star_h11)
            
#             if not np.any(np.isnan(star_h11)):
            fesc_m_bin.append(np.divide(flux_h11,star_h11))

#             plt.plot(_x,np.divide(flux_h11,star_h11), 'o', c=color[i], ls=ls[i])
# _________________________________________________________________________________________________
#         y=np.mean(flux_m_bin, axis=0)/np.mean(star_m_bin, axis=0)
#         y=np.mean( np.divide(flux_m_bin,star_m_bin), axis=0)
#         y=np.nanmean(fesc_m_bin, axis=0)        
#         plt.ylabel("R200 flux over Stellar flux [unknow unit]")
# _________________________________________________________________________________________________
        y=np.nanmean(flux_m_bin, axis=0)
        plt.ylabel("R200 flux [unknow unit]")

# _________________________________________________________________________________________________
#         y=np.nanmean(star_m_bin, axis=0)
#         plt.ylabel("Stellar flux [phot/s]")
# _________________________________________________________________________________________________        
        ls=["-",":"]
        print("%s -> %.1e <Mh< %.1e"%(ls[i], Mbins[i],Mbins[i+1]))
        print("%s -> %s"%(color[istep], labels[istep]))
        plt.plot(_x,y, 'o', ls=ls[i],  label=lab, c=color[istep])


plt.yscale("log")
plt.xlabel("r/R200")
plt.xlim(0,1)
# plt.ylim(1e-57,1e-51

