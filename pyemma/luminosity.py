
# coding: utf-8

# In[ ]:

import os, sys
import numpy as np
from scipy import integrate


# In[5]:

def mag2flux(mag):
    flux = np.power(10,-mag/2.5)
    return flux


# In[ ]:

def flux2mag(flux):
    mag = -2.5*np.log10(flux)
    return mag


# In[5]:

def readStarburst(convert=False):

    filename="/home/deparis/jupyter/Starburst99/fig1e.dat.txt"
    #reading header
    with open(filename, 'r') as file:
        file.readline()
        file.readline()
        name=file.readline().split()[2:]

    age = [	float(n[:-3])*1e6 for n in name]

    #loading data
    data = np.loadtxt(filename,unpack=True, skiprows=3)
    data[1:] = np.power(10,data[1:])

    if convert :
        print ("converting unit to SI")
        data[0]  *= 1e-10 #A->m
        data[1:] *= 1e-7/1e-10 #log(Erg/A) -> J/m/s

    return age,data[0], data[1:]


# In[10]:

def getM1600(_x0, spectremodeleenergparsparAngstrom):
    """
    Pour calculer M1600 pour un spectre modèle défini par spectremodeleenergparsparAngstrom 
    et ses songeurs d’onde _x0:

    flux de ref F0: 3631 Jy
    Jyhzcm=1.e-23; //  erg s-1 Hz-1 cm-2
    Jyhzm=1.e-19; // erg s-1 Hz-1 m-2
    // so for absolute mags the 0 point is Jyhzm * 4pir^2 where r = 10 pc
    AB0pointabs=3631.*Jyhzm*4.*pi*(10.*ParSec)^2; // erg/s/Hz // cause ParSec
    """


    F0 = 3631 #Jy
    Jyhzcm=1.e-23; #  erg s-1 Hz-1 cm-2
    Jyhzm=1.e-19; # erg s-1 Hz-1 m-2
    #so for absolute mags the 0 point is Jyhzm * 4pir^2 where r = 10 pc
    parsec=3.08567758e16
    AB0pointabs=F0*Jyhzm*4.*np.pi*(10.*parsec)**2; # erg/s/Hz // cause ParSec

    """
    // or flam = Fnu * (-c/lam^2)
    // so lets make a reference spectrum corresponding to the 0-pointABabsloulte
    AB0pointref=AB0pointabs*(c*1.e10/((_x0)^2)); // in erg/s/A (thats why we have to put c in Angstrom/s)
    _x0 c’est les longueurs d’onde en Angstrom de ton spectre de pop stellaire
    """
    c=299792458 #light speed
    AB0pointref=AB0pointabs*(c*1.e10/((_x0)**2)); # in erg/s/A (thats why we have to put c in Angstrom/s)

    """
    fi1600=_x0*0.;   // creation d’un filtre
    fi1600(180:205)=1.;   // le filtre vaut 1 entre 1500 et 1600 Angstrom
    MAB1600=-2.5*log10(integ(fi1600*spectremodeleenergparsparAngstrom,_x0)/integ(fi1600*AB0pointref,_x0));
    MAB1600;

    A répéter pour tous les spectres modèle. 
    """

    mask = np.where( (_x0>=1500) & (_x0<=1600))#le filtre vaut 1 entre 1500 et 1600 Angstrom
    MAB1600=-2.5*np.log10(integrate.trapz(spectremodeleenergparsparAngstrom[mask],
                                          _x0[mask])/integrate.trapz(AB0pointref[mask],
                                        _x0[mask]))

    return MAB1600


# In[1]:

def getModel():
    age, wavelength, spectre = readStarburst()

    n = len(age)
    modelMag = np.zeros(n)

    for i in range(n):
        modelMag[i] = getM1600(wavelength, spectre[i])

    return  modelMag, age


# In[9]:

def sumMag(mass, age, modelmag, modelage):
    mags=np.interp(age,modelage,modelmag)
    res=np.sum(mass*np.power(10,-mags/2.5))
    return -2.5*np.log10(res)


# In[14]:

def get_all_flux_1600(stars,current_time,unit_mass):
    """
    get luminous flux of stars in M1600 band
    """
    
    modelmag, modelage = getModel()
    age = current_time - stars.age.data
    mass = stars.mass.data*unit_mass/1.9891e30/1e6
    mags=np.interp(age,modelage,modelmag)
    flux = mass * mag2flux(mags)
    
    stars.flux_1600=flux


# In[6]:

def get_all_flux_UV(stars,current_time,unit_mass):
    """
    get luminous flux of stars in UV
    """
    
    def getflux_UV(age):
        tlife = 3.673e6 #yr
        E0 = 3.399e16 #phot/s/kg
        y=np.ones(len(age)) *E0
        y[age>tlife] *= np.power(age[age>tlife]/tlife ,-4.)
        return y
        
    mass = stars.mass.data*unit_mass # kg
    flux = getflux_UV(current_time - stars.age.data) # #phot/s/kg
    stars.flux_UV = mass * flux #phot/s


# In[3]:

def get_tot_egy(age, tlife):
    """
    return the total energy emmited by a source 
    """
    def f(t, tlife):
        if t<=tlife :
            return 1.
        else:
            return np.power(t/tlife ,-4.)
    
    return integrate.quad(f,0,age, args=(tlife))[0]

