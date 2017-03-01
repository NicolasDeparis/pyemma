
# coding: utf-8

# In[1]:

get_ipython().magic('matplotlib notebook')
import matplotlib.pyplot as plt
import numpy as np
import sys

from scipy import integrate
from scipy import spatial

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


# # READ

# In[3]:

run1=io.Run(runset.runs[1].folder)


# In[4]:

run1.step_00015.grid.x.read()
run1.step_00015.grid.y.read()
run1.step_00015.grid.z.read()
run1.step_00015.grid.l.read()
run1.step_00015.grid.field_d.read()
# run1.step_00015.grid.z_last_xion.read()
# run1.step_00015.grid.xion.read()


# In[6]:

# run1.step_00015.star.x.read()
# run1.step_00015.star.y.read()
# run1.step_00015.star.z.read()
run1.step_00015.star.mass.read()
run1.step_00015.star.age.read()


# In[5]:

run2=io.Run(runset.runs[0].folder)


# In[329]:

run2.step_00017.grid.x.read()
run2.step_00017.grid.y.read()
run2.step_00017.grid.z.read()
run2.step_00017.grid.l.read()
run2.step_00017.grid.field_d.read()
run2.step_00017.grid.z_last_xion.read()
run2.step_00017.grid.xion.read()


# In[7]:

run2.step_00017.star.x.read()
run2.step_00017.star.y.read()
run2.step_00017.star.z.read()
run2.step_00017.star.mass.read()
run2.step_00017.star.age.read()


# In[ ]:

x=cat1.x
y=cat1.y
z=cat1.z
cat1_tree = spatial.KDTree( np.transpose( [x,y,z] ))


# In[46]:

run3=io.Run(runset.runs[4].folder)


# In[88]:

run3.step_00019.grid.x.read()
run3.step_00019.grid.y.read()
run3.step_00019.grid.z.read()
run3.step_00019.grid.l.read()
run3.step_00019.grid.field_d.read()
run3.step_00019.grid.z_last_xion.read()


# In[47]:

run3.step_00019.star.x.read()
run3.step_00019.star.y.read()
run3.step_00019.star.z.read()
run3.step_00019.star.mass.read()
run3.step_00019.star.age.read()


# In[8]:

run4=io.Run(runset.runs[2].folder)


# In[51]:

run4.step_00012.grid.x.read()
run4.step_00012.grid.y.read()
run4.step_00012.grid.z.read()
run4.step_00012.grid.l.read()
run4.step_00012.grid.field_d.read()
run4.step_00012.grid.z_last_xion.read()


# In[10]:

run4.step_00012.star.x.read()
run4.step_00012.star.y.read()
run4.step_00012.star.z.read()
run4.step_00012.star.mass.read()
run4.step_00012.star.age.read()


# In[319]:

run5=io.Run(runset.runs[8].folder)


# In[320]:

run5.step_00015.grid.x.read()
run5.step_00015.grid.y.read()
run5.step_00015.grid.z.read()
run5.step_00015.grid.l.read()
run5.step_00015.grid.field_d.read()
run5.step_00015.grid.z_last_xion.read()


# In[321]:

run5.step_00015.star.x.read()
run5.step_00015.star.y.read()
run5.step_00015.star.z.read()
run5.step_00015.star.mass.read()
run5.step_00015.star.age.read()


# # FOF

# In[5]:

cat1=fof.Fof(runset.runs[1].folder+"data/",15,8)
cat1.getNfofTot()
cat1.get_masst()

ob=run1.param.info.ob
om=run1.param.info.om
cat1.get_R200(ob, om, 8)
cat1.get_stars(run1.step_00015.star)
cat1.get_cells(run1.step_00015.grid)


# In[ ]:

age=run1.step_00015.star.age.data
mass=run1.step_00015.star.mass.data
tlife = run1.param.run.tlife
cat1.get_integ_egy(age,mass,tlife, force=1)


# In[54]:

cat2=fof.Fof(runset.runs[0].folder+"data/",17,8)

# cat2.write_fofin()
# cat2.write_infosim(256,8)
# cat2.gen(8)

cat2.getNfofTot()
cat2.get_masst()

ob=run2.param.info.ob
om=run2.param.info.om
cat2.get_R200(ob, om, 8)
cat2.get_stars(run2.step_00017.star)
cat2.get_cells(run2.step_00017.grid)


# In[55]:

cat3=fof.Fof(runset.runs[4].folder+"data/",19,8)

# cat3.write_fofin()
# cat3.write_infosim(256,8)
# cat3.gen(8)

cat3.getNfofTot()
cat3.get_masst()

ob=run3.param.info.ob
om=run3.param.info.om
cat3.get_R200(ob, om, 8)
cat3.get_stars(run3.step_00019.star)
cat3.get_cells(run3.step_00019.grid)


# In[56]:

cat4=fof.Fof(runset.runs[2].folder+"data/",12,8)

# cat4.write_fofin()
# cat4.write_infosim(256,8)
# cat4.gen(8)

cat4.getNfofTot()
cat4.get_masst()

ob=run4.param.info.ob
om=run4.param.info.om
cat4.get_R200(ob, om, 8)
cat4.get_stars(run4.step_00012.star)
cat4.get_cells(run4.step_00012.grid)


# In[322]:

cat5=fof.Fof(runset.runs[8].folder+"data/",15,8)

cat5.write_fofin()
cat5.write_infosim(256,8)
cat5.gen(8)

cat5.getNfofTot()
cat5.get_masst()

ob=run5.param.info.ob
om=run5.param.info.om
cat5.get_R200(ob, om, 8)
cat5.get_stars(run5.step_00015.star)
cat5.get_cells(run5.step_00015.grid)


# # PLOT

# In[177]:

n=3

halo_num1=np.argsort(cat1.npart)[-n]
stars1=cat1.stars[halo_num1]

halo_num2=np.argsort(cat2.npart)[-n]
stars2=cat2.stars[halo_num2]

halo_num3=np.argsort(cat3.npart)[-n]
stars3=cat3.stars[halo_num3]

halo_num4=np.argsort(cat4.npart)[-n]
stars4=cat4.stars[halo_num4]

fig=plt.figure()
ax = fig.add_subplot(1,1,1)
ax.add_patch(plt.Circle((cat1.x[halo_num1],cat1.y[halo_num1]),cat1.R200[halo_num1],fill=False))
ax.add_patch(plt.Circle((cat2.x[halo_num2],cat2.y[halo_num2]),cat2.R200[halo_num2],fill=False))
ax.add_patch(plt.Circle((cat3.x[halo_num3],cat3.y[halo_num3]),cat3.R200[halo_num3],fill=False))
ax.add_patch(plt.Circle((cat4.x[halo_num4],cat4.y[halo_num4]),cat4.R200[halo_num4],fill=False))

# x=run4.step_00012.star.x.data[stars4]
# y=run4.step_00012.star.y.data[stars4]
# plt.plot(x,y,'b.')

x=run3.step_00019.star.x.data[stars3]
y=run3.step_00019.star.y.data[stars3]
plt.plot(x,y,'y.')

x=run2.step_00017.star.x.data[stars2]
y=run2.step_00017.star.y.data[stars2]
plt.plot(x,y,'r.')

x=run1.step_00015.star.x.data[stars1]
y=run1.step_00015.star.y.data[stars1]
plt.plot(x,y,'b.')


# # MOMENT D'INERTIE

# In[192]:

xc=cat1.x[halo_num1]
yc=cat1.y[halo_num1]
zc=cat1.z[halo_num1]

x=run1.step_00015.star.x.data[stars1]
y=run1.step_00015.star.y.data[stars1]
z=run1.step_00015.star.z.data[stars1]
m=run1.step_00015.star.mass.data[stars1]

r= np.sqrt( pow(x-xc,2)+pow(y-yc,2)+pow(y-yc,2))
print(np.sum(r*m ) )

xc=cat2.x[halo_num2]
yc=cat2.y[halo_num2]
zc=cat2.z[halo_num2]

x=run2.step_00017.star.x.data[stars2]
y=run2.step_00017.star.y.data[stars2]
z=run2.step_00017.star.z.data[stars2]
m=run2.step_00017.star.mass.data[stars2]

r= np.sqrt( pow(x-xc,2)+pow(y-yc,2)+pow(y-yc,2))
print( np.sum(r*m ) )

xc=cat3.x[halo_num3]
yc=cat3.y[halo_num3]
zc=cat3.z[halo_num3]

x=run3.step_00019.star.x.data[stars3]
y=run3.step_00019.star.y.data[stars3]
z=run3.step_00019.star.z.data[stars3]
m=run3.step_00019.star.mass.data[stars3]

r= np.sqrt( pow(x-xc,2)+pow(y-yc,2)+pow(y-yc,2))
print( np.sum(r*m ) )


# xc=cat4.x[halo_num4]
# yc=cat4.y[halo_num4]
# zc=cat4.z[halo_num4]

# x=run4.step_00012.star.x.data[stars4]
# y=run4.step_00012.star.y.data[stars4]
# z=run4.step_00012.star.z.data[stars4]
# m=run4.step_00012.star.mass.data[stars4]

# r= np.sqrt( pow(x-xc,2)+pow(y-yc,2)+pow(y-yc,2))
# print( np.sum(r*m ) )


# # STAR DENSITY PROFIL

# In[185]:

n=3
nbins=32
plt.figure()


################################################################################
# xc=cat2.x[halo_num2]
# yc=cat2.y[halo_num2]
# zc=cat2.z[halo_num2]
# halo_num1=cat1_tree.query((xc,yc,zc))[1]

arg= np.argsort(cat3.npart)
halo_num3=arg[-n]

xc=cat3.x[halo_num3]
yc=cat3.y[halo_num3]
zc=cat3.z[halo_num3]
print(xc,yc,zc)

stars3=cat3.stars[halo_num3]

x=run3.step_00019.star.x.data[stars3]
y=run3.step_00019.star.y.data[stars3]
z=run3.step_00019.star.z.data[stars3]
m=run3.step_00019.star.mass.data[stars3]
print(np.sum(m))

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/512")

################################################################################

arg= np.argsort(cat2.npart)
halo_num2=arg[-n]

stars2=cat2.stars[halo_num2]

xc=cat2.x[halo_num2]
yc=cat2.y[halo_num2]
zc=cat2.z[halo_num2]
print(xc,yc,zc)

x=run2.step_00017.star.x.data[stars2]
y=run2.step_00017.star.y.data[stars2]
z=run2.step_00017.star.z.data[stars2]
m=run2.step_00017.star.mass.data[stars2]
print(np.sum(m))

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *3.

plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/64")

################################################################################
# xc=cat2.x[halo_num2]
# yc=cat2.y[halo_num2]
# zc=cat2.z[halo_num2]
# halo_num1=cat1_tree.query((xc,yc,zc))[1]


arg= np.argsort(cat1.npart)
halo_num1=arg[-n]

xc=cat1.x[halo_num1]
yc=cat1.y[halo_num1]
zc=cat1.z[halo_num1]
print(xc,yc,zc)

stars1=cat1.stars[halo_num1]

x=run1.step_00015.star.x.data[stars1]
y=run1.step_00015.star.y.data[stars1]
z=run1.step_00015.star.z.data[stars1]
m=run1.step_00015.star.mass.data[stars1]
print(np.sum(m))

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/8")

################################################################################



plt.legend()
plt.yscale("log")

plt.xlabel("radius")
plt.ylabel("Stellar density")


# # LUMINOSITY PROFIL

# In[ ]:

run1.param.info.om


# In[15]:

def a2t_quad(az):

    o_m=0.316
    o_v=1.-o_m

    H0=67
    H0 = H0*1e3/1e6/3.08567758e16*(365*24*3600) # Hubble constant in SI unit

    def f(a):
        return a/ np.sqrt(o_m*a + o_v*a**4 )

    return 1./H0 * integrate.quad(f,0,az)[0]


# In[323]:

cur_step=run1.step_00015
# cur_step=run2.step_00017
# cur_step=run3.step_00019
# cur_step=run4.step_00012
# cur_step=run5.step_00015

stars=cur_step.star
a=cur_step.star.x._tsim
print(a)
t=a2t_quad(0.15)
unit_m=run2.param.info.unit_mass
luminosity.get_all_flux(stars,t,a,unit_m)


# In[7]:

n=4
nbins=32

arg= np.argsort(cat1.npart)
halo_num1=arg[-n]

xc=cat1.x[halo_num1]
yc=cat1.y[halo_num1]
zc=cat1.z[halo_num1]

stars1=cat1.stars[halo_num1]

x=run1.step_00015.star.x.data[stars1]
y=run1.step_00015.star.y.data[stars1]
z=run1.step_00015.star.z.data[stars1]


r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

flux = run1.step_00015.star.flux[stars1]

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=flux)
h12,_=np.histogram(r,bins=bins,weights=flux*flux)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *0


h11=luminosity.flux2mag(h11)
err1=luminosity.flux2mag(err1)

plt.figure()
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/8")

print(luminosity.flux2mag(np.sum(flux)))

########################################################################################

cur_step=run2.step_00017
arg= np.argsort(cat2.npart)
halo_num2=arg[-n]

xc=cat2.x[halo_num2]
yc=cat2.y[halo_num2]
zc=cat2.z[halo_num2]

stars2=cat2.stars[halo_num2]

x=cur_step.star.x.data[stars2]
y=cur_step.star.y.data[stars2]
z=cur_step.star.z.data[stars2]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

# bins=np.linspace(np.min(r),np.max(r),nbins+1)
# _x=(bins[1:]+bins[:-1])/2

flux =cur_step.star.flux[stars2]

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=flux)
h12,_=np.histogram(r,bins=bins,weights=flux*flux)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *0

h11=luminosity.flux2mag(h11)
err1=luminosity.flux2mag(err1)
print(luminosity.flux2mag(np.sum(flux)))

plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/64")
########################################################################################

cur_step=run3.step_00019
arg= np.argsort(cat3.npart)
halo_num3=arg[-n]

xc=cat3.x[halo_num3]
yc=cat3.y[halo_num3]
zc=cat3.z[halo_num3]

stars3=cat3.stars[halo_num3]

x=cur_step.star.x.data[stars3]
y=cur_step.star.y.data[stars3]
z=cur_step.star.z.data[stars3]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

# bins=np.linspace(np.min(r),np.max(r),nbins+1)
# _x=(bins[1:]+bins[:-1])/2

flux =cur_step.star.flux[stars3]

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=flux)
h12,_=np.histogram(r,bins=bins,weights=flux*flux)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *0

h11=luminosity.flux2mag(h11)
err1=luminosity.flux2mag(err1)
print(luminosity.flux2mag(np.sum(flux)))

plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/512")



plt.gca().invert_yaxis()
# plt.yscale("log")
plt.xlabel("radius [box unit]")
plt.ylabel("Mag")
plt.legend()


# # AGE PROFIL

# In[186]:

n=3
nbins=32
plt.figure()


################################################################################

arg= np.argsort(cat3.npart)
halo_num3=arg[-n]

stars3=cat3.stars[halo_num3]

xc=cat3.x[halo_num3]
yc=cat3.y[halo_num3]
zc=cat3.z[halo_num3]
print(xc,yc,zc)

x=run3.step_00019.star.x.data[stars3]
y=run3.step_00019.star.y.data[stars3]
z=run3.step_00019.star.z.data[stars3]
m=run3.step_00019.star.age.data[stars3]
print(np.sum(m))

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *3.

h11=h11/n1
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/64")

################################################################################

arg= np.argsort(cat2.npart)
halo_num2=arg[-n]

stars2=cat2.stars[halo_num2]

xc=cat2.x[halo_num2]
yc=cat2.y[halo_num2]
zc=cat2.z[halo_num2]
print(xc,yc,zc)

x=run2.step_00017.star.x.data[stars2]
y=run2.step_00017.star.y.data[stars2]
z=run2.step_00017.star.z.data[stars2]
m=run2.step_00017.star.age.data[stars2]
print(np.sum(m))

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *3.

h11=h11/n1
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/64")

################################################################################
# xc=cat2.x[halo_num2]
# yc=cat2.y[halo_num2]
# zc=cat2.z[halo_num2]
# halo_num1=cat1_tree.query((xc,yc,zc))[1]

arg= np.argsort(cat1.npart)
halo_num1=arg[-n]

xc=cat1.x[halo_num1]
yc=cat1.y[halo_num1]
zc=cat1.z[halo_num1]
print(xc,yc,zc)

stars1=cat1.stars[halo_num1]

x=run1.step_00015.star.x.data[stars1]
y=run1.step_00015.star.y.data[stars1]
z=run1.step_00015.star.z.data[stars1]
m=run1.step_00015.star.age.data[stars1]
print(np.sum(m))

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

n1,_=np.histogram(r,bins=bins)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1)  /np.sqrt(n1) *3.

h11=h11/n1
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/8")

plt.legend()
plt.yscale("log")


# # DENSITY PROFIL

# In[172]:

n=3
nbins=64
plt.figure()
################################################################################

arg= np.argsort(cat1.npart)
halo_num1=arg[-n]


xc=cat1.x[halo_num1]
yc=cat1.y[halo_num1]
zc=cat1.z[halo_num1]

cells=cat1.cells[halo_num1]

cur_step=run1.step_00015
x=cur_step.grid.x.data[cells]
y=cur_step.grid.y.data[cells]
z=cur_step.grid.z.data[cells]
l=cur_step.grid.l.data[cells]
d=cur_step.grid.field_d.data[cells]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))
dv=np.power(0.5,3*l)
m=d*dv

bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2

n1,_=np.histogram(r,bins=bins,weights=dv)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *3.

h11=h11/n1
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/8")

################################################################################

arg= np.argsort(cat2.npart)
halo_num2=arg[-n]


xc=cat2.x[halo_num2]
yc=cat2.y[halo_num2]
zc=cat2.z[halo_num2]

cells=cat2.cells[halo_num2]

x=run2.step_00017.grid.x.data[cells]
y=run2.step_00017.grid.y.data[cells]
z=run2.step_00017.grid.z.data[cells]
l=run2.step_00017.grid.l.data[cells]
d=run2.step_00017.grid.field_d.data[cells]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))
dv=np.power(0.5,3*l)
m=d*dv

n1,_=np.histogram(r,bins=bins,weights=dv)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *3.

h11=h11/n1
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/64")

################################################################################

arg= np.argsort(cat3.npart)
halo_num3=arg[-n]

xc=cat3.x[halo_num3]
yc=cat3.y[halo_num3]
zc=cat3.z[halo_num3]

cells=cat3.cells[halo_num3]

x=run3.step_00019.grid.x.data[cells]
y=run3.step_00019.grid.y.data[cells]
z=run3.step_00019.grid.z.data[cells]
l=run3.step_00019.grid.l.data[cells]
d=run3.step_00019.grid.field_d.data[cells]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))
dv=np.power(0.5,3*l)
m=d*dv

n1,_=np.histogram(r,bins=bins,weights=dv)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *3.

h11=h11/n1
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/512")



plt.legend()
plt.xlabel("Radius [Box unit]")
plt.ylabel("Overdensity")
plt.yscale("log")


# # Z_REIO PROFIL

# In[187]:

n=3
nbins=32
plt.figure()
################################################################################

arg= np.argsort(cat1.npart)
halo_num1=arg[-n]

xc=cat1.x[halo_num1]
yc=cat1.y[halo_num1]
zc=cat1.z[halo_num1]

cells=cat1.cells[halo_num1]

cur_step=run1.step_00015
x=cur_step.grid.x.data[cells]
y=cur_step.grid.y.data[cells]
z=cur_step.grid.z.data[cells]
l=cur_step.grid.l.data[cells]
zxion=cur_step.grid.z_last_xion.data[cells]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))
dv=np.power(0.5,3*l)
m=zxion*dv


bins=np.linspace(np.min(r),np.max(r),nbins+1)
_x=(bins[1:]+bins[:-1])/2


n1,_=np.histogram(r,bins=bins,weights=dv)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *3.

h11=h11/n1
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/8")

################################################################################

arg= np.argsort(cat2.npart)
halo_num2=arg[-n]


xc=cat2.x[halo_num2]
yc=cat2.y[halo_num2]
zc=cat2.z[halo_num2]

cells=cat2.cells[halo_num2]

x=run2.step_00017.grid.x.data[cells]
y=run2.step_00017.grid.y.data[cells]
z=run2.step_00017.grid.z.data[cells]
l=run2.step_00017.grid.l.data[cells]
zxion=run2.step_00017.grid.z_last_xion.data[cells]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))
dv=np.power(0.5,3*l)
m=zxion*dv

n1,_=np.histogram(r,bins=bins,weights=dv)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *3.

h11=h11/n1
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/64")

################################################################################

arg= np.argsort(cat3.npart)
halo_num3=arg[-n]

xc=cat3.x[halo_num3]
yc=cat3.y[halo_num3]
zc=cat3.z[halo_num3]

cells=cat3.cells[halo_num3]

x=run3.step_00019.grid.x.data[cells]
y=run3.step_00019.grid.y.data[cells]
z=run3.step_00019.grid.z.data[cells]
l=run3.step_00019.grid.l.data[cells]
zxion=run3.step_00019.grid.z_last_xion.data[cells]

r= np.sqrt( np.power(x-xc,2)+np.power(y-yc,2)+np.power(z-zc,2))
dv=np.power(0.5,3*l)
m=zxion*dv

n1,_=np.histogram(r,bins=bins,weights=dv)
h11,_=np.histogram(r,bins=bins,weights=m)
h12,_=np.histogram(r,bins=bins,weights=m*m)
err1 = np.sqrt(h12/n1 - h11*h11/n1/n1) /np.sqrt(n1) *3.

h11=h11/n1
plt.errorbar(_x[1:],h11[1:],yerr=err1[1:],label="1/512")



plt.legend()
plt.xlabel("Radius [Box unit]")
plt.ylabel("Reionisation redshift")
# plt.xscale("log")


# # LUMINOSITY FUNCTION

# In[2]:

flux = run1.step_00015.star.flux

n=cat1.nfoftot
flux_tot1=np.zeros(n)
for i in range(n):
    stars=cat1.stars[i]
    flux_tot1[i]=np.sum(flux[stars])    
mag1=luminosity.flux2mag(flux_tot1[flux_tot1!=0])


# In[51]:

flux = run2.step_00017.star.flux

n=cat2.nfoftot
flux_tot2=np.zeros(n)
for i in range(n):
    stars=cat2.stars[i]
    flux_tot2[i]=np.sum(flux[stars])
mag2=luminosity.flux2mag(flux_tot2[flux_tot2!=0])


# In[77]:

flux = run3.step_00019.star.flux

n=cat3.nfoftot
flux_tot3=np.zeros(n)
for i in range(n):
    stars=cat3.stars[i]
    flux_tot3[i]=np.sum(flux[stars])
mag3=luminosity.flux2mag(flux_tot3[flux_tot3!=0])


# In[78]:

flux = run4.step_00012.star.flux

n=cat4.nfoftot
flux_tot4=np.zeros(n)
for i in range(n):
    stars=cat4.stars[i]
    flux_tot4[i]=np.sum(flux[stars])
mag4=luminosity.flux2mag(flux_tot4[flux_tot4!=0])


# In[324]:

flux = run5.step_00015.star.flux

n=cat5.nfoftot
flux_tot5=np.zeros(n)
for i in range(n):
    stars=cat5.stars[i]
    flux_tot5[i]=np.sum(flux[stars])
mag5=luminosity.flux2mag(flux_tot5[flux_tot5!=0])


# In[50]:

nbins=16
plt.figure()


mag_min= min(np.min(mag2),np.min(mag4))
mag_max= max(np.max(mag2),np.max(mag4))

bins=np.linspace(mag_min,mag_max,nbins+1)
_x=(bins[1:]+bins[:-1])/2
dx=np.diff(_x)

################################################################################

n1,_=np.histogram(mag1,bins=bins)
err1=np.sqrt(n1)
box_size=run1.param.info.box_size_hm1_Mpc
plt.errorbar(_x[1:],n1[1:]/dx/ (box_size/0.67)**3,label="1/8")

################################################################################

n1,_=np.histogram(mag2,bins=bins)

box_size=run2.param.info.box_size_hm1_Mpc
err1=np.sqrt(n1)
plt.errorbar(_x[1:],n1[1:]/dx/ (box_size/0.67)**3,label="1/64")

################################################################################

n1,_=np.histogram(mag3,bins=bins)

box_size=run3.param.info.box_size_hm1_Mpc
err1=np.sqrt(n1)
plt.errorbar(_x[1:],n1[1:]/dx/ (box_size/0.67)**3,label="1/512")

################################################################################

# n1,_=np.histogram(mag4,bins=bins)

# box_size=run4.param.info.box_size_hm1_Mpc
# err1=np.sqrt(n1)
# plt.errorbar(_x[1:],n1[1:]/dx/ (box_size/0.67)**3,label="1/64 SN off")

################################################################################

# n1,_=np.histogram(mag5,bins=bins)

# box_size=run5.param.info.box_size_hm1_Mpc
# err1=np.sqrt(n1)
# plt.errorbar(_x[1:],n1[1:]/dx/ (box_size/0.67)**3,label="1/64 SN off")

################################################################################

x,y=observations.luminosity_function_fit(6)
plt.plot(x,y,'k--')

plt.ylim(1e-4,1e-0)
plt.yscale("log", nonposy='mask')
plt.legend(loc=0)


# # LUMINOSITY FUNCTION OF STELLAR MASS

# In[194]:

plt.figure()

########################################################################################################

a=run1.step_00015.star.x._tsim
print(a)
a=0.15
t=a2t_quad(a)

mass = run1.step_00015.star.mass.data
age = t - run1.step_00015.star.age.data

n=cat1.nfoftot
mass_tot=np.zeros(n)
age_tot=np.zeros(n)
for i in range(n):
    stars=cat1.stars[i]
    
    mass_tot[i]=np.sum(mass[stars])
    age_tot[i]=np.mean(age[stars])

mass_tot = mass_tot/1.9891e30*run1.param.info.unit_mass
# mag1=luminosity.flux2mag(flux_tot1[flux_tot1!=0])

# plt.plot(age_tot[flux_tot1!=0],mag1,'.',label="1/8")
# plt.plot(age_tot,mass_tot,'.')
plt.plot(mass_tot[flux_tot1!=0],mag1,'.',label="1/8")

########################################################################################################


cur_step=run2.step_00017
a=cur_step.star.mass._tsim
print(a)
t=a2t_quad(a)

mass = cur_step.star.mass.data
age = t - cur_step.star.age.data

n=cat2.nfoftot
mass_tot=np.zeros(n)
age_tot=np.zeros(n)
for i in range(n):
    stars=cat2.stars[i]
    
    mass_tot[i]=np.sum(mass[stars])
    age_tot[i]=np.mean(age[stars])

    
mass_tot = mass_tot/1.9891e30*run2.param.info.unit_mass

# plt.plot(age_tot[flux_tot2!=0],mag2,'.',label='1./64')
# plt.plot(age_tot,mass_tot,'.')

plt.plot(mass_tot[flux_tot2!=0],mag2,'.',label='1/64')

########################################################################################################

cur_step=run3.step_00019
a=cur_step.star.mass._tsim
# A=0.15
print(a)
t=a2t_quad(a)

mass =cur_step.star.mass.data
age = t - cur_step.star.age.data


n=cat3.nfoftot
mass_tot=np.zeros(n)
age_tot=np.zeros(n)
for i in range(n):
    stars=cat3.stars[i]
    
    mass_tot[i]=np.sum(mass[stars])
    age_tot[i]=np.mean(age[stars])

    
mass_tot = mass_tot/1.9891e30*run3.param.info.unit_mass

# plt.plot(age_tot[flux_tot2!=0],mag2,'.',label='1./64')
# plt.plot(age_tot,mass_tot,'.')

plt.plot(mass_tot[flux_tot3!=0],mag3,'.',label='1/512')

########################################################################################################

# cur_step=run4.step_00012
# a=cur_step.star.mass._tsim
# # A=0.15
# print(a)
# t=a2t_quad(a)

# mass =cur_step.star.mass.data
# age = t - cur_step.star.age.data


# n=cat4.nfoftot
# mass_tot=np.zeros(n)
# age_tot=np.zeros(n)
# for i in range(n):
#     stars=cat4.stars[i]
    
#     mass_tot[i]=np.sum(mass[stars])
#     age_tot[i]=np.mean(age[stars])

    
# mass_tot = mass_tot/1.9891e30*run4.param.info.unit_mass

# plt.plot(age_tot[flux_tot2!=0],mag2,'.',label='1./64')
# plt.plot(age_tot,mass_tot,'.')

# plt.plot(mass_tot[flux_tot4!=0],mag4,'.',label='1/64 SN off')

########################################################################################################
plt.legend(loc=0)
plt.xscale('log')
plt.gca().invert_yaxis()
# plt.yscale('log')
plt.xlabel("Stellar mass")
plt.ylabel("Mag")


# # STELLAR MASS FUNCTION OF HALO MASS

# In[195]:

level=8
nptot=2**(3*level)
ob=run2.param.info.ob
om=run2.param.info.om
part_mass=(1.-ob/om)/nptot
print(part_mass/1.9891e30*run2.param.info.unit_mass)

plt.figure()

##########################################################################""
mass=cat1.npart*part_mass                
mass = mass/1.9891e30*run1.param.info.unit_mass

n=cat1.nfoftot
mass_tot=np.zeros(n)
for i in range(n):
    stars=cat1.stars[i]    
    mass_tot[i]=np.sum(run1.step_00015.star.mass.data[stars])
mass_tot = mass_tot/1.9891e30*run1.param.info.unit_mass
    
plt.plot(mass,mass_tot,'.')

##########################################################################""

mass=cat2.npart*part_mass                
mass = mass/1.9891e30*run2.param.info.unit_mass

n=cat2.nfoftot
mass_tot=np.zeros(n)
for i in range(n):
    stars=cat2.stars[i]    
    mass_tot[i]=np.sum(run2.step_00017.star.mass.data[stars])
mass_tot = mass_tot/1.9891e30*run2.param.info.unit_mass
    
plt.plot(mass,mass_tot,'.')

##########################################################################""

mass=cat3.npart*part_mass                
mass = mass/1.9891e30*run3.param.info.unit_mass

n=cat3.nfoftot
mass_tot=np.zeros(n)
for i in range(n):
    stars=cat3.stars[i]    
    mass_tot[i]=np.sum(run3.step_00019.star.mass.data[stars])
mass_tot = mass_tot/1.9891e30*run3.param.info.unit_mass
    
plt.plot(mass,mass_tot,'.')

##########################################################################""

# mass=cat4.npart*part_mass                
# mass = mass/1.9891e30*run4.param.info.unit_mass

# n=cat4.nfoftot
# mass_tot=np.zeros(n)
# for i in range(n):
#     stars=cat4.stars[i]    
#     mass_tot[i]=np.sum(run4.step_00012.star.mass.data[stars])
# mass_tot = mass_tot/1.9891e30*run4.param.info.unit_mass
    
# plt.plot(mass,mass_tot,'.')
##########################################################################""
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Halo mass [M0]')
plt.ylabel('Stellar mass [M0]')


# In[10]:

plt.figure()
plt.plot(cat1.npart, cat1.integ_egy, '.')
plt.xscale("log")
plt.yscale("log")


# # Z_REIO FUNCTION OF MAG

# In[307]:

level=8
nptot=2**(3*level)
ob=run2.param.info.ob
om=run2.param.info.om
part_mass=(1.-ob/om)/nptot
# print(part_mass/1.9891e30*run2.param.info.unit_mass)

plt.figure()

##########################################################################""

mass=cat1.npart*part_mass                
mass = mass/1.9891e30*run1.param.info.unit_mass

    
mean_zxion=np.zeros(cat1.nfoftot)
v=np.zeros(cat1.nfoftot)

for halo_num in range(cat1.nfoftot):
    
    cells=cat1.cells[halo_num]
    l=run1.step_00015.grid.l.data[cells]
    x=run1.step_00015.grid.z_last_xion.data[cells]

    v[halo_num] = np.sum( np.power(0.5,3*l) )
    mean_zxion[halo_num] = np.sum( x* np.power(0.5,3*l) )
# plt.plot(mag1[mag1>-10],(mean_zxion[flux_tot1!=0]/v[flux_tot1!=0])[mag1>-10],'o',label="1/8")
plt.plot(mag1,(mean_zxion[flux_tot1!=0]/v[flux_tot1!=0]),'.',label="1/8")


# mass=cat1.npart*part_mass                
# mass = mass/1.9891e30*run4.param.info.unit_mass
# plt.plot(mass[flux_tot1!=0],mean_zxion[flux_tot1!=0]/v[flux_tot1!=0],'o',label="1/8")

##########################################################################""

mass=cat2.npart*part_mass                
mass = mass/1.9891e30*run2.param.info.unit_mass

    
mean_zxion=np.zeros(cat2.nfoftot)
v=np.zeros(cat2.nfoftot)

for halo_num in range(cat2.nfoftot):
    
    cells=cat2.cells[halo_num]
    l=run2.step_00017.grid.l.data[cells]
    x=run2.step_00017.grid.z_last_xion.data[cells]

    v[halo_num] = np.sum( np.power(0.5,3*l) )
    mean_zxion[halo_num] = np.sum( x* np.power(0.5,3*l) )
    
plt.plot(mag2,mean_zxion[flux_tot2!=0]/v[flux_tot2!=0],'r.', label="1/64")



# mass=cat2.npart*part_mass                
# mass = mass/1.9891e30*run4.param.info.unit_mass
# plt.plot(mass[flux_tot2!=0],mean_zxion[flux_tot2!=0]/v[flux_tot2!=0],'r.', label="1/64")


##########################################################################""

mass=cat3.npart*part_mass                
mass = mass/1.9891e30*run3.param.info.unit_mass

    
mean_zxion=np.zeros(cat3.nfoftot)
v=np.zeros(cat3.nfoftot)

for halo_num in range(cat3.nfoftot):
    
    cells=cat3.cells[halo_num]
    l=run3.step_00019.grid.l.data[cells]
    x=run3.step_00019.grid.z_last_xion.data[cells]

    v[halo_num] = np.sum( np.power(0.5,3*l) )
    mean_zxion[halo_num] = np.sum( x* np.power(0.5,3*l) )
    
plt.plot(mag3,mean_zxion[flux_tot3!=0]/v[flux_tot3!=0],'.', label="1/512")

##########################################################################""
# plt.xscale('log')
# plt.yscale('log')
plt.legend()
plt.gca().invert_xaxis()
# plt.xlabel('Halo mass [M0]')
plt.xlabel('Mag')
plt.ylabel('Mean z of reionization')


# In[446]:

level=8
nptot=2**(3*level)
ob=run2.param.info.ob
om=run2.param.info.om
part_mass=(1.-ob/om)/nptot
# print(part_mass/1.9891e30*run2.param.info.unit_mass)

plt.figure()

##########################################################################""

mass=cat1.npart*part_mass                
mass = mass/1.9891e30*run1.param.info.unit_mass
    
mean_zxion=np.zeros(cat1.nfoftot)
v=np.zeros(cat1.nfoftot)

for halo_num in range(cat1.nfoftot):
    
    cells=cat1.cells[halo_num]
    l=run1.step_00015.grid.l.data[cells]
    x=run1.step_00015.grid.z_last_xion.data[cells]

    v[halo_num] = np.sum( np.power(0.5,3*l) )
    mean_zxion[halo_num] = np.sum( x* np.power(0.5,3*l) )

plt.plot(mass[flux_tot1!=0],mean_zxion[flux_tot1!=0]/v[flux_tot1!=0],'bo',label="1/8")

##########################################################################""

mass=cat2.npart*part_mass                
mass = mass/1.9891e30*run2.param.info.unit_mass
    
mean_zxion=np.zeros(cat2.nfoftot)
v=np.zeros(cat2.nfoftot)

for halo_num in range(cat2.nfoftot):
    
    cells=cat2.cells[halo_num]
    l=run2.step_00017.grid.l.data[cells]
    x=run2.step_00017.grid.z_last_xion.data[cells]

    v[halo_num] = np.sum( np.power(0.5,3*l) )
    mean_zxion[halo_num] = np.sum( x* np.power(0.5,3*l) )
    
plt.plot(mass[flux_tot2!=0],mean_zxion[flux_tot2!=0]/v[flux_tot2!=0],'r.', label="1/64")


##########################################################################""

mass=cat3.npart*part_mass                
mass = mass/1.9891e30*run3.param.info.unit_mass

mean_zxion=np.zeros(cat3.nfoftot)
v=np.zeros(cat3.nfoftot)

for halo_num in range(cat3.nfoftot):
    
    cells=cat3.cells[halo_num]
    l=run3.step_00019.grid.l.data[cells]
    x=run3.step_00019.grid.z_last_xion.data[cells]

    v[halo_num] = np.sum( np.power(0.5,3*l) )
    mean_zxion[halo_num] = np.sum( x* np.power(0.5,3*l) )
    
plt.plot(mass[flux_tot3!=0],mean_zxion[flux_tot3!=0]/v[flux_tot3!=0],'y.', label="1/512")
##########################################################################""
plt.xscale('log')
# plt.yscale('log')
plt.legend()
# plt.xlabel('Halo mass [M0]')
plt.xlabel('Mass [Mo]')
plt.ylabel('Mean z of reionization')


# # TOTAL FLUX FUNCTION OF MASS

# In[49]:

plt.figure()

level=32
nptot=2**(3*level)
ob=run1.param.info.ob
om=run1.param.info.om
part_mass=(1.-ob/om)/nptot

mass=cat1.npart*part_mass                
mass = mass/1.9891e30*run1.param.info.unit_mass

mass_min= np.min(mass)
mass_max= np.max(mass)

nbins=16
# bins=np.linspace(mass_min,mass_max,nbins+1)
bins=np.logspace(np.log10(mass_min),np.log10(mass_max),nbins+1)
_x=(bins[1:]+bins[:-1])/2
dx=np.diff(_x)

n0,_=np.histogram(mass,bins=bins)
n1,_=np.histogram(mass,bins=bins, weights=flux_tot1)
n2,_=np.histogram(mass,bins=bins, weights=flux_tot1*flux_tot1)
err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3. *0

n1=n1/n0
n1=luminosity.flux2mag(n1)
err=luminosity.flux2mag(err)

box_size=run1.param.info.box_size_hm1_Mpc
plt.errorbar(_x[1:],n1[1:],yerr=err[1:],label="1/8")

################################################################################

mass=cat2.npart*part_mass                
mass = mass/1.9891e30*run2.param.info.unit_mass

n0,_=np.histogram(mass,bins=bins)
n1,_=np.histogram(mass,bins=bins, weights=flux_tot2)
n2,_=np.histogram(mass,bins=bins, weights=flux_tot2*flux_tot2)
err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3. *0

n1=n1/n0
n1=luminosity.flux2mag(n1)
err=luminosity.flux2mag(err)

box_size=run1.param.info.box_size_hm1_Mpc
plt.errorbar(_x[1:],n1[1:],yerr=err[1:],label="1/64")

################################################################################

mass=cat3.npart*part_mass                
mass = mass/1.9891e30*run3.param.info.unit_mass

n0,_=np.histogram(mass,bins=bins)
n1,_=np.histogram(mass,bins=bins, weights=flux_tot3)
n2,_=np.histogram(mass,bins=bins, weights=flux_tot3*flux_tot3)
err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3. *0

n1=n1/n0
n1=luminosity.flux2mag(n1)
err=luminosity.flux2mag(err)

box_size=run1.param.info.box_size_hm1_Mpc
plt.errorbar(_x[1:],n1[1:],yerr=err[1:],label="1/512")


plt.legend(loc=0)
plt.xscale('log')
# plt.yscale('log',nonposy='mask')
plt.xlabel("Halo Mass [M0]")
plt.ylabel("Mag")
plt.gca().invert_yaxis()


# In[448]:

plt.figure()

level=8
nptot=2**(3*level)
ob=run1.param.info.ob
om=run1.param.info.om
part_mass=(1.-ob/om)/nptot

################################################################################

mass=cat3.npart*part_mass                
mass = mass/1.9891e30*run3.param.info.unit_mass

mag=luminosity.flux2mag(flux_tot3)
plt.plot(mass,mag,'.',label="1/512")

################################################################################

mass=cat2.npart*part_mass                
mass = mass/1.9891e30*run2.param.info.unit_mass

mag=luminosity.flux2mag(flux_tot2)
plt.plot(mass,mag,'.',label="1/64")

################################################################################

mass=cat1.npart*part_mass                
mass = mass/1.9891e30*run1.param.info.unit_mass

mag=luminosity.flux2mag(flux_tot1)
plt.plot(mass,mag,'.',label="1/8")

################################################################################


plt.legend(loc=0)
plt.xscale('log')
plt.xlabel("Halo Mass [M0]")
plt.ylabel("Mag")
plt.gca().invert_yaxis()

plt.axvline(2e8,  color='k', ls='--')
plt.axhline(-11.6,  color='k', ls='--')


# In[ ]:

level=8
nptot=2**(3*level)
ob=run1.param.info.ob
om=run1.param.info.om
part_mass=(1.-ob/om)/nptot

mass=cat1.npart*part_mass                
mass = mass/1.9891e30*run1.param.info.unit_mass

mag=luminosity.flux2mag(flux_tot1)

mask= (mass<8e7) & (mag<-13.7)
print(cat1.x[mask])
print(cat1.y[mask])
print(cat1.z[mask])


# BARYON

# In[205]:

baryon_mass1=np.zeros(cat1.nfoftot)

for halo_num in range(cat1.nfoftot):
    
    cells=cat1.cells[halo_num]
    l=run1.step_00015.grid.l.data[cells]
    d=run1.step_00015.grid.field_d.data[cells]

    baryon_mass1[halo_num] = np.sum( d* np.power(0.5,3*l) )
        
baryon_mass1 = baryon_mass1/1.9891e30*run1.param.info.unit_mass


# In[27]:

baryon_mass2=np.zeros(cat2.nfoftot)

for halo_num in range(cat2.nfoftot):
    
    cells=cat2.cells[halo_num]
    l=run2.step_00017.grid.l.data[cells]
    d=run2.step_00017.grid.field_d.data[cells]

    baryon_mass2[halo_num] = np.sum( d* np.power(0.5,3*l) )
baryon_mass2 = baryon_mass2/1.9891e30*run2.param.info.unit_mass


# In[213]:

baryon_mass3=np.zeros(cat3.nfoftot)

for halo_num in range(cat3.nfoftot):
    
    cells=cat3.cells[halo_num]
    l=run3.step_00019.grid.l.data[cells]
    d=run3.step_00019.grid.field_d.data[cells]

    baryon_mass3[halo_num] = np.sum( d* np.power(0.5,3*l) )
baryon_mass3 = baryon_mass3/1.9891e30*run3.param.info.unit_mass


# In[28]:

baryon_mass4=np.zeros(cat4.nfoftot)

for halo_num in range(cat4.nfoftot):
    
    cells=cat4.cells[halo_num]
    l=run4.step_00012.grid.l.data[cells]
    d=run4.step_00012.grid.field_d.data[cells]

    baryon_mass4[halo_num] = np.sum( d* np.power(0.5,3*l) )
baryon_mass4 = baryon_mass4/1.9891e30*run4.param.info.unit_mass


# In[ ]:

plt.figure()
mass_min=np.min(baryon_mass1[baryon_mass1!=0])
mass_max=np.max(baryon_mass1)

nbins=32
bins=np.logspace(np.log10(mass_min),np.log10(mass_max),nbins+1)
_x=(bins[1:]+bins[:-1])/2


n0,_=np.histogram(baryon_mass1, bins=bins)
plt.plot(_x,n0)

n0,_=np.histogram(baryon_mass2, bins=bins)
plt.plot(_x,n0)

n0,_=np.histogram(baryon_mass3, bins=bins)
plt.plot(_x,n0)

n0,_=np.histogram(baryon_mass4, bins=bins)
plt.plot(_x,n0)

plt.xscale('log')
plt.yscale('log')


# In[210]:

level=8
nptot=2**(3*level)
ob=run1.param.info.ob
om=run1.param.info.om
part_mass=(1.-ob/om)/nptot

mass_tot1=cat1.npart*part_mass
mass_tot1 = mass_tot1/1.9891e30*run1.param.info.unit_mass


# In[30]:

level=8
nptot=2**(3*level)
ob=run2.param.info.ob
om=run2.param.info.om
part_mass=(1.-ob/om)/nptot

mass_tot2=cat2.npart*part_mass
mass_tot2=mass_tot2/1.9891e30*run2.param.info.unit_mass


# In[211]:

level=8
nptot=2**(3*level)
ob=run3.param.info.ob
om=run3.param.info.om
part_mass=(1.-ob/om)/nptot

mass_tot3=cat3.npart*part_mass
mass_tot3=mass_tot3/1.9891e30*run3.param.info.unit_mass


# In[31]:

level=8
nptot=2**(3*level)
ob=run4.param.info.ob
om=run4.param.info.om
part_mass=(1.-ob/om)/nptot

mass_tot4=cat4.npart*part_mass
mass_tot4=mass_tot4/1.9891e30*run4.param.info.unit_mass


# In[207]:

n=cat1.nfoftot
stars_mass_tot1=np.zeros(n)
for i in range(n):
    stars=cat1.stars[i]    
    stars_mass_tot1[i]=np.sum(run1.step_00015.star.mass.data[stars])
stars_mass_tot1 = stars_mass_tot1/1.9891e30*run1.param.info.unit_mass


# In[32]:

n=cat2.nfoftot
stars_mass_tot2=np.zeros(n)
for i in range(n):
    stars=cat2.stars[i]    
    stars_mass_tot2[i]=np.sum(run2.step_00017.star.mass.data[stars])
stars_mass_tot2 = stars_mass_tot2/1.9891e30*run2.param.info.unit_mass


# In[215]:

n=cat3.nfoftot
stars_mass_tot3=np.zeros(n)
for i in range(n):
    stars=cat3.stars[i]    
    stars_mass_tot3[i]=np.sum(run3.step_00019.star.mass.data[stars])
stars_mass_tot3 = stars_mass_tot3/1.9891e30*run3.param.info.unit_mass
print(stars_mass_tot3)


# In[33]:

n=cat4.nfoftot
stars_mass_tot4=np.zeros(n)
for i in range(n):
    stars=cat4.stars[i]    
    stars_mass_tot4[i]=np.sum(run4.step_00012.star.mass.data[stars])
stars_mass_tot4 = stars_mass_tot4/1.9891e30*run4.param.info.unit_mass


# In[ ]:

mag1=luminosity.flux2mag(flux_tot1)
mag2=luminosity.flux2mag(flux_tot2)


# # LUMINOSITY FUNCTION OF BARYONIC FRACTION

# In[217]:

plt.figure()

baryon_frac1= (baryon_mass1+stars_mass_tot1)/mass_tot1
baryon_frac2= (baryon_mass2+stars_mass_tot2)/mass_tot2
baryon_frac3= (baryon_mass3+stars_mass_tot3)/mass_tot3
# baryon_frac4= (baryon_mass4+stars_mass_tot4)/mass_tot4


nbins=32
mass_min= np.nanmin(baryon_frac2[flux_tot2!=0])
mass_max= np.nanmax(baryon_frac2[flux_tot2!=0])

bins=np.logspace(np.log10(mass_min),np.log10(mass_max),nbins+1)
dx=np.diff(bins)
_x=(bins[1:]+bins[:-1])/2


#########################################################################################
mask=flux_tot1!=0
w=flux_tot1[mask]
n0,_=np.histogram(baryon_frac1[mask],bins=bins)
n1,_=np.histogram(baryon_frac1[mask],bins=bins, weights=w)
n2,_=np.histogram(baryon_frac1[mask],bins=bins, weights=w*w)

err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3.*0

_y=luminosity.flux2mag(n1/n0)
err=luminosity.flux2mag(err)

# plt.plot(_x,_y,'bo--', label='1/8')
plt.errorbar(_x,_y,yerr=err,label="1/8")

#########################################################################################

mask=flux_tot2!=0
w=flux_tot2[mask]
n0,_=np.histogram(baryon_frac2[mask],bins=bins)
n1,_=np.histogram(baryon_frac2[mask],bins=bins, weights=w)
n2,_=np.histogram(baryon_frac2[mask],bins=bins, weights=w*w)

err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3.*0

_y=luminosity.flux2mag(n1/n0)
err=luminosity.flux2mag(err)

# plt.plot(_x,_y,'bo--', label='1/8')
plt.errorbar(_x,_y,yerr=err,label="1/64")

#########################################################################################

mask=flux_tot3!=0
w=flux_tot3[mask]
n0,_=np.histogram(baryon_frac3[mask],bins=bins)
n1,_=np.histogram(baryon_frac3[mask],bins=bins, weights=w)
n2,_=np.histogram(baryon_frac3[mask],bins=bins, weights=w*w)

err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3.*0

_y=luminosity.flux2mag(n1/n0)
err=luminosity.flux2mag(err)

# plt.plot(_x,_y,'bo--', label='1/8')
plt.errorbar(_x,_y,yerr=err,label="1/512")

#########################################################################################

# mask=flux_tot4!=0
# w=flux_tot4[mask]
# n0,_=np.histogram(baryon_frac4[mask],bins=bins)
# n1,_=np.histogram(baryon_frac4[mask],bins=bins, weights=w)
# n2,_=np.histogram(baryon_frac4[mask],bins=bins, weights=w*w)

# err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3.*0

# _y=luminosity.flux2mag(n1/n0)
# err=luminosity.flux2mag(err)

# # plt.plot(_x,_y,'bo--', label='1/8')
# plt.errorbar(_x,_y,yerr=err,label="1/64 SN off")

#########################################################################################
# plt.plot(baryon_frac2,mag2,'ro')

plt.ylabel("Mag")
plt.xlabel("Baryonic fraction")
plt.xscale('log')
plt.gca().invert_yaxis()
plt.legend(loc=0)

# universal baryonic fraction
ob=run2.param.info.ob
om=run2.param.info.om
plt.axvline(ob/om,  color='k', ls='--')


# In[218]:

plt.figure()

#########################################################################################

# baryon_frac4= (baryon_mass4+stars_mass_tot4)/mass_tot4
# mask=flux_tot4!=0
# _x=baryon_frac4[mask]
# _y=luminosity.flux2mag(flux_tot4[mask])

# plt.plot(_x,_y,'y.', label='1/64 SN off')

#########################################################################################

baryon_frac3= (baryon_mass3+stars_mass_tot3)/mass_tot3
mask=flux_tot3!=0
_x=baryon_frac3[mask]
_y=luminosity.flux2mag(flux_tot3[mask])

plt.plot(_x,_y,'y.', label='1/512')

#########################################################################################

baryon_frac2= (baryon_mass2+stars_mass_tot2)/mass_tot2
mask=flux_tot2!=0
_x=baryon_frac2[mask]
_y=luminosity.flux2mag(flux_tot2[mask])

plt.plot(_x,_y,'r.', label='1/64')
#########################################################################################

baryon_frac1= (baryon_mass1+stars_mass_tot1)/mass_tot1
mask=flux_tot1!=0
_x=baryon_frac1[mask]
_y=luminosity.flux2mag(flux_tot1[mask])

plt.plot(_x,_y,'bo', label='1/8')

#########################################################################################

# universal baryonic fraction
ob=run2.param.info.ob
om=run2.param.info.om
plt.axvline(ob/om,  color='k', ls='--')

#########################################################################################

plt.ylabel("Mag")
plt.xlabel("Baryonic fraction")
plt.xscale('log')
plt.gca().invert_yaxis()
plt.legend(loc=2)


# # BARYNIC FRACTION FUNCTION OF HALO MASS

# In[303]:

plt.figure()

nbins=32

mask=flux_tot1!=0
mask = np.ones(mass_tot1.shape, dtype=np.bool)

mass_min= np.nanmin(mass_tot1[mask])
mass_max= np.nanmax(mass_tot1[mask])

bins=np.logspace(np.log10(mass_min),np.log10(mass_max),nbins+1)

dx=np.diff(bins)
_x=(bins[1:]+bins[:-1])/2

#########################################################################################

baryon_frac1= (baryon_mass1+stars_mass_tot1)/mass_tot1

w=baryon_frac1[mask]

n0,_=np.histogram(mass_tot1[mask],bins=bins)
n1,_=np.histogram(mass_tot1[mask],bins=bins, weights=w)
n2,_=np.histogram(mass_tot1[mask],bins=bins, weights=w*w)
err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3

_y=n1/n0

plt.errorbar(_x,_y,yerr=err,label="1/8", fmt='.-')

#########################################################################################

baryon_frac2= (baryon_mass2+stars_mass_tot2)/mass_tot2

w=baryon_frac2[mask]

n0,_=np.histogram(mass_tot2[mask],bins=bins)
n1,_=np.histogram(mass_tot2[mask],bins=bins, weights=w)
n2,_=np.histogram(mass_tot2[mask],bins=bins, weights=w*w)
err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3

_y=n1/n0

plt.errorbar(_x,_y,yerr=err,label="1/8", fmt='.-')

#########################################################################################

baryon_frac3= (baryon_mass3+stars_mass_tot3)/mass_tot3

w=baryon_frac3[mask]

n0,_=np.histogram(mass_tot3[mask],bins=bins)
n1,_=np.histogram(mass_tot3[mask],bins=bins, weights=w)
n2,_=np.histogram(mass_tot3[mask],bins=bins, weights=w*w)
err = np.sqrt(n2/n0 - n1*n1/n0/n0) /np.sqrt(n0) *3

_y=n1/n0

plt.errorbar(_x,_y,yerr=err,label="1/8", fmt='.-')

#########################################################################################

# universal baryonic fraction
ob=run1.param.info.ob
om=run1.param.info.om
plt.axhline(ob/om,  color='k', ls='--')

plt.xlabel("Halo mass")
plt.ylabel("Baryon fraction")
plt.xscale('log')
plt.yscale('log')


# In[220]:

plt.figure()

plt.plot(mass_tot3,baryon_mass3/mass_tot3,'.')

plt.plot(mass_tot2,baryon_mass2/mass_tot2,'.')
plt.plot(mass_tot1,baryon_mass1/mass_tot1,'.')

# plt.plot(mass_tot1,baryon_frac1,'.')
# plt.plot(mass_tot2,baryon_frac2,'.')
# plt.plot(mass_tot3,baryon_frac3,'.')
# plt.plot(mass_tot4,baryon_frac4,'.')

# universal baryonic fraction
ob=run2.param.info.ob
om=run2.param.info.om
plt.axhline(ob/om,  color='k', ls='--')

plt.xlabel("Halo mass")
plt.ylabel("Baryon fraction")
plt.xscale('log')
plt.yscale('log')


# # HALO MAP

# In[440]:

fig=plt.figure()
# ax = fig.add_subplot(1,1,1)
# for i in range(cat1.nfoftot):
#     ax.add_patch(plt.Circle((cat1.x[i],cat1.y[i]),cat1.R200[i],fill=False))

# mask1=[flux_tot1!=0][mag1>-10]

plt.plot(cat1.x,cat1.y,'k.')
plt.plot(cat1.x[flux_tot1!=0],cat1.y[flux_tot1!=0],'bo')

x= np.copy(cat1.x)
y= np.copy(cat1.y)

tmp = np.zeros(cat1.nfoftot)
mask = (flux_tot1!=0) & (mass_tot1<1e8)
x[ np.logical_not(mask)  ]=0
y[ np.logical_not(mask)  ]=0

mask = mag1<-11.6
x[ np.logical_not(mask)  ]=0
y[ np.logical_not(mask)  ]=0

plt.plot(x,y,'ro')


# # IONIZATION STATE INDICATOR

# $log_{10}\left(\left(1+\delta\right) \left( \frac{{x_{HII}}^2}{x_{HI}} \right) \right)$ function of $\left(1+\delta\right) $

# In[345]:

plt.figure
_x=run1.step_00015.grid.field_d.data
_y=run1.step_00015.grid.xion.data
_l=run1.step_00015.grid.l.data
dv=np.power(0.5,3*_l)

plt.figure()

mask = np.where(_y!=1.)
_x = _x[mask]
_y = _y[mask]
_y = _x * np.power(_y,2)/(1.-_y)

dv = dv[mask]
V=np.sum(dv)

n=256

x=np.log10(_x)
y=np.log10(_y)

H,yedges,xedges=np.histogram2d(y,x,bins=n, weights=dv/V)
X, Y = np.meshgrid(xedges, yedges)
H = np.ma.masked_where(H==0,H) 

plt.pcolormesh(X, Y, np.log(H), cmap='hot')

plt.xlabel(r'$log10 \left( 1+\delta \right)$')
plt.ylabel(r'$log10\left(\left(1+\delta\right) \left( \frac{{x_{HII}}^2}{x_{HI}} \right) \right)$')

cb=plt.colorbar()
cb.set_label('log10 of number of occurence')

plt.xlim(xedges[0], xedges[-1])
plt.ylim(yedges[0], yedges[-1])


# In[346]:

plt.figure
_x=run2.step_00017.grid.field_d.data
_y=run2.step_00017.grid.xion.data
_l=run2.step_00017.grid.l.data
dv=np.power(0.5,3*_l)

plt.figure()

mask = np.where(_y!=1.)
_x = _x[mask]
_y = _y[mask]
_y = _x * np.power(_y,2)/(1.-_y)

dv = dv[mask]
V=np.sum(dv)

n=256

x=np.log10(_x)
y=np.log10(_y)

H,yedges,xedges=np.histogram2d(y,x,bins=n, weights=dv/V)
X, Y = np.meshgrid(xedges, yedges)
H = np.ma.masked_where(H==0,H) 

plt.pcolormesh(X, Y, np.log(H), cmap='hot')

plt.xlabel(r'$log10 \left( 1+\delta \right)$')
plt.ylabel(r'$log10\left(\left(1+\delta\right) \left( \frac{{x_{HII}}^2}{x_{HI}} \right) \right)$')

cb=plt.colorbar()
cb.set_label('log10 of number of occurence')

plt.xlim(xedges[0], xedges[-1])
plt.ylim(yedges[0], yedges[-1])


# In[354]:

khi1=np.zeros(cat1.nfoftot)

cur_step= run1.step_00015
for halo_num in range(cat1.nfoftot):
    
    cells=cat1.cells[halo_num]
    l=cur_step.grid.l.data[cells]
    d=cur_step.grid.field_d.data[cells]
    x=cur_step.grid.xion.data[cells]
    
    dv=np.power(0.5,3*l)    
    V = np.sum(dv)
    
    khi_loc = d * np.power(x,2)/(1.-x)
    khi1[halo_num] = np.sum( khi_loc*dv ) /V
    


# In[368]:

khi2=np.zeros(cat2.nfoftot)

cur_step= run2.step_00017
for halo_num in range(cat2.nfoftot):
    
    cells=cat2.cells[halo_num]
    l=cur_step.grid.l.data[cells]
    d=cur_step.grid.field_d.data[cells]
    x=cur_step.grid.xion.data[cells]
    
    dv=np.power(0.5,3*l)    
    V = np.sum(dv)
    
    khi_loc = d * np.power(x,2)/(1.-x)
    khi2[halo_num] = np.sum( khi_loc*dv ) /V


# In[367]:

khi3=np.zeros(cat3.nfoftot)

cur_step= run3.step_00019
for halo_num in range(cat3.nfoftot):
    
    cells=cat3.cells[halo_num]
    l=cur_step.grid.l.data[cells]
    d=cur_step.grid.field_d.data[cells]
    x=cur_step.grid.xion.data[cells]
    
    dv=np.power(0.5,3*l)    
    V = np.sum(dv)
    
    khi_loc = d * np.power(x,2)/(1.-x)
    khi3[halo_num] = np.sum( khi_loc*dv ) /V


# In[375]:

plt.figure()


plt.plot(mass_tot1,khi1,'.')
plt.plot(mass_tot2,khi2,'.')
# plt.plot(mass_tot3,khi3,'.')
plt.xscale('log')

# plt.plot(mag1,khi1[flux_tot1!=0],'.')
# plt.plot(mag2,khi2[flux_tot2!=0],'.')
# plt.plot(mag3,khi3[flux_tot3!=0],'.')
# plt.gca().invert_xaxis()


plt.xlabel("Halo mass")
plt.ylabel("Mean ionisation state")
plt.yscale('log')


# # X

# In[390]:

wv1=np.zeros(cat1.nfoftot)
wm1=np.zeros(cat1.nfoftot)

cur_step= run1.step_00015
for halo_num in range(cat1.nfoftot):
    
    cells=cat1.cells[halo_num]
    l=cur_step.grid.l.data[cells]
    d=cur_step.grid.field_d.data[cells]
    x=1-cur_step.grid.xion.data[cells]
    
    dv=np.power(0.5,3*l)    
    V = np.sum(dv)
    M = np.sum(d*dv)
    
    wv1[halo_num]=np.sum(x*dv)/V
    wm1[halo_num]=np.sum(x*d*dv)/M


# In[391]:

cur_step= run2.step_00017
cur_cat = cat2


wv2=np.zeros(cur_cat.nfoftot)
wm2=np.zeros(cur_cat.nfoftot)

for halo_num in range(cur_cat.nfoftot):
    
    cells=cur_cat.cells[halo_num]
    l=cur_step.grid.l.data[cells]
    d=cur_step.grid.field_d.data[cells]
    x=1-cur_step.grid.xion.data[cells]
    
    dv=np.power(0.5,3*l)    
    V = np.sum(dv)
    M = np.sum(d*dv)
    
    wv2[halo_num]=np.sum(x*dv)/V
    wm2[halo_num]=np.sum(x*d*dv)/M


# In[392]:

cur_step= run3.step_00019
cur_cat = cat3

wv3=np.zeros(cur_cat.nfoftot)
wm3=np.zeros(cur_cat.nfoftot)

for halo_num in range(cur_cat.nfoftot):
    
    cells=cur_cat.cells[halo_num]
    l=cur_step.grid.l.data[cells]
    d=cur_step.grid.field_d.data[cells]
    x=1-cur_step.grid.xion.data[cells]
    
    dv=np.power(0.5,3*l)    
    V = np.sum(dv)
    M = np.sum(d*dv)
    
    wv3[halo_num]=np.sum(x*dv)/V
    wm3[halo_num]=np.sum(x*d*dv)/M


# In[396]:

plt.figure()

plt.plot(wv1,wm1,'.')
plt.plot(wv2,wm2,'.')
plt.plot(wv3,wm3,'.')

plt.xlabel("xion volume weighted")
plt.ylabel("xion mass weighted")
plt.xscale('log')
plt.yscale('log')


# In[398]:

plt.figure()


plt.plot(mass_tot1,wm1/wv1,'.')
plt.plot(mass_tot2,wm2/wv2,'.')
plt.plot(mass_tot3,wm3/wv3,'.')

plt.xlabel("Halo mass [Mo]")
plt.ylabel("xion mass weighted over xion volume weighted")
plt.xscale('log')
plt.yscale('log')


# In[467]:

plt.figure()

plt.plot(stars_mass_tot1,wm1/wv1,'.')
plt.plot(stars_mass_tot2,wm2/wv2,'.')
plt.plot(stars_mass_tot3,wm3/wv3,'.')

plt.xlabel("Stellar mass [Mo]")
plt.ylabel("xion mass weighted over xion volume weighted")
plt.xscale('log')
plt.yscale('log')


# # SFR

# In[32]:

def t2a(t):
    """
        convert time to expansion factor
    """
    n = 100
    A = np.arange(n) / float(n)
    
    T=np.zeros(n)
    for i in range(n):
        T[i] = a2t_quad(A[i])
        
    return np.interp(t,T,A)


def getFromPop(stars,param,mask=None, Nbins=16):

    if mask==None:
        mask=np.ones(stars.age.data.shape[0], dtype=np.bool)
    
    z=np.zeros(Nbins)
    sfr=np.zeros(Nbins)
        
    ages =  stars.age.data[mask]
    
    age_min = np.min(ages)
    age_max = np.max(ages)
    bins_edges = np.linspace(age_min,age_max,Nbins+1)
    
    cur_t=a2t_quad(stars.age._tsim)
    masses= stars.mass.data[mask]  
    masses[cur_t- ages>param.run.tlife_SN] /= (1-param.run.ejecta_proportion)    
    masses = masses/1.9891e30*param.info.unit_mass

    sfr,bin_edges=np.histogram(ages, bins=bins_edges, weights=masses)
    sfr/=np.diff(bin_edges)

    a=t2a(bin_edges)
    z_all=1./a-1
    z=(z_all[1:]+z_all[:-1])/2

    return z,sfr


# In[31]:

run1.step_00015.star.age._tsim=0.15
z, sfr = getFromPop(run1.step_00015.star, run1.param)
sfr/=np.power(run1.param.info.box_size_hm1_Mpc/ run1.param.info.H0*100,3)


plt.figure()
plt.plot(z,sfr, label=lab)
    
observations.sfr()
plt.yscale('log')


# # SFR FUNCTION OF HALO MASS

# In[1]:

level=8
nptot=2**(3*level)
ob=run1.param.info.ob
om=run1.param.info.om
part_mass=(1.-ob/om)/nptot

mass_tot1=cat1.npart*part_mass
mass_tot1=mass_tot1/1.9891e30*run1.param.info.unit_mass

mass_min=np.min(mass_tot1)
mass_max=np.max(mass_tot1)
nbins = 1
bins=np.logspace(np.log10(mass_min),np.log10(mass_max),nbins+1)
_x=(bins[1:]+bins[:-1])/2

stars=np.zeros(nbins,dtype=np.object)
for i in range(nbins):    
    mask=  (mass_tot1>=bins[i]) & (mass_tot1<bins[i+1])            
    star = np.concatenate(cat1.stars[mask]).astype(np.int)  
    stars[i]=starCX_FREEZE
    
plt.figure()    
for i in range(nbins):
    
    z, sfr = getFromPop(run1.step_00015.star, run1.param, stars[i])    
    sfr/=np.power(run1.param.info.box_size_hm1_Mpc/ run1.param.info.H0*100,3)
    lab = "%.1e < M < %.1e"%(bins[i],bins[i+1])
    plt.plot(z,sfr, label=lab)
    
observations.sfr()
plt.yscale('log')
plt.legend()


# In[ ]:

cat1.get


# In[682]:

level=8
nptot=2**(3*level)
ob=run2.param.info.ob
om=run2.param.info.om
part_mass=(1.-ob/om)/nptot

mass_tot1=cat1.npart*part_mass
mass_tot1=mass_tot1/1.9891e30*run1.param.info.unit_mass

mass_min=np.min(mass_tot1)
mass_max=np.max(mass_tot1)
nbins = 8
bins=np.logspace(np.log10(mass_min),np.log10(mass_max),nbins+1)
_x=(bins[1:]+bins[:-1])/2


stars=np.zeros(nbins,dtype=np.object)
for i in range(nbins):    
    mask=  (mass_tot1>=bins[i]) & (mass_tot1<bins[i+1])            
    star = np.concatenate(cat1.stars[mask]).astype(np.int)
    stars[i]=star

plt.figure()

_y=np.zeros(nbins)

for i in range(nbins):
    y=0
    print(stars[i].shape)
    for j in range(stars[i].shape[0]):
        y+=luminosity.get_tot_egy(run1.step_00015.star.age.data[stars[i][j]  ] , run1.param)
    _y=y
    
plt.plot(_x, _y, 'o-')
    

    
    



# In[73]:

n=1
halo_num1=np.argsort(cat1.npart)[-n]

#center
xc = cat1.x[halo_num1]
yc = cat1.y[halo_num1]
zc = cat1.z[halo_num1]

#halo cells
cells=cat1.cells[halo_num1]

grid_x=run1.step_00015.grid.x.data[cells]
grid_y=run1.step_00015.grid.y.data[cells]
grid_z=run1.step_00015.grid.z.data[cells]
grid_l=run1.step_00015.grid.l.data[cells]

#radius
R = np.sqrt(  np.power(grid_x-xc,2)
            + np.power(grid_y-yc,2)
            + np.power(grid_z-zc,2) )

R200=cat1.R200[halo_num1]
dr=np.power(0.5,np.min(grid_l))

#selection
mask= (R>R200-dr) & (R<=R200)

x=grid_x[mask]
y=grid_y[mask]
z=grid_z[mask]


# In[9]:

data = np.loadtxt("../jupyter/healpix")
n_healpix=data[0]
x_healpix=data[1:n_healpix+1]
y_healpix=data[n_healpix+1:2*n_healpix+1]
z_healpix=data[2*n_healpix+1:3*n_healpix+1]


# In[44]:

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(9,9))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_healpix*R200, y_healpix*R200, z_healpix*R200, '.')


# In[47]:

halo_num = halo_num1

cells = cat1.cells[halo_num]
R200=cat1.R200[halo_num]

cur_step = run1.step_00015
x_grid = cur_step.grid.x.data[cells]
y_grid = cur_step.grid.y.data[cells]
z_grid = cur_step.grid.z.data[cells]
l = cur_step.grid.l.data[cells]

dx = np.power(0.5,l) 
idx = np.zeros(n, dtype=np.int32)
n=x_healpix.shape[0]
for i in range(n):
    x=(x_healpix *R200/2 + cat1.x[halo_num])[i]
    y=(y_healpix *R200/2 + cat1.y[halo_num])[i]
    z=(z_healpix *R200/2 + cat1.z[halo_num])[i]

    xm = x >= x_grid
    xp = x <  x_grid+dx
    ym = y >= y_grid
    yp = y <  y_grid+dx
    zm = z >= z_grid
    zp = z <  z_grid+dx
    
    res = xm & xp & ym & yp & zm & zp
    idx[i] = np.int32(np.where(res==True))

d = cur_step.grid.field_d.data[cells]
print(np.mean(d[idx]))


# In[22]:




# In[164]:

from scipy import interpolate

grid=np.array( (x_grid,y_grid,z_grid) ).T
val = run1.step_00015.grid.field_d.data[cells]
points = (x_healpix *R200 + cat1.x[halo_num1], y_healpix*R200 + cat1.y[halo_num1], z_healpix*R200 + cat1.y[halo_num1]) 
interpolate.griddata(grid, val, points, method='nearest')


# In[48]:

get_ipython().magic('pinfo np.argsort')

