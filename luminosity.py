# coding: utf-8

import os, sys
import numpy as np
from scipy import integrate

def mag2flux(mag):
    flux = np.power(10,-mag/2.5)
    return flux

def flux2mag(flux):
    mag = -2.5*np.log10(flux)
    return mag

def readStarburst_old(convert=False, filename=''):

    if( filename=='' ):
        filename="/astro/home/nicolas.gillet/WORK/analyse_simu/fig7e.dat"
       
    print( 'Reading %s'%filename )
    #reading header
    with open(filename, 'r') as file:
        file.readline()
        file.readline()
        name=file.readline().split()[2:]

    age = np.array( [float(n[:-3])*1e6 for n in name] )

    #loading data
    data = np.loadtxt( filename, unpack=True, skiprows=3 )
    data[1:] = np.power(10,data[1:])

    if convert :
        print ("converting unit to SI")
        data[0]  *= 1e-10 #A->m
        data[1:] *= 1e-7/1e-10 #log(Erg/A) -> J/m/s

    return age,data[0], data[1:]

def readStarburst(filename=''):

    # file = "../jupyter/Starburst99/Salpeter/SALP_AUBERT.spectrum1"
    # file = "../jupyter/Starburst99/Topheavy/TH_AUBERT.spectrum1"
    # file = "../jupyter/Starburst99/Kroupa/KROUPA_DOM.spectrum1"
    if( filename=='' ):
        filename = "/astro/home/nicolas.gillet/WORK/analyse_simu/TH_AUBERT.spectrum1"
        
    print( 'Reading %s'%filename )

    data = np.loadtxt(filename,unpack=True,skiprows=6)
    time = np.unique(data[0,:])

    spectrum = []
    for cur_time in time:
            mask=data[0,:]==cur_time
            cur_data = data[:,mask]
            wavelenght = cur_data[1]
            spectrum.append(np.power(10,cur_data[2]))

    return time, wavelenght, np.array(spectrum)

def readBPASS(filename=''):
    """
    read BPASS spectrum
    time [y]
    wavelenght [A]
    spectrum [erg.s-1.A-1]
    """
    
    if filename=='':
        filename = "/astro/home/nicolas.gillet/WORK/analyse_simu/BPASSv2_imf135_300/OUTPUT_POP/spectra-bin.z001.dat"
        
    print( 'Reading %s'%filename )
    
    data = np.loadtxt(filename, unpack=True)
    
    spectrum = data[1:,:] * 3.8274e26 * 1.e7 ### convert for Luminosity_sun.A-1 to erg.s-1.A-1
    wavelenght = data[0]
    time = 10**( 6 +0.1*np.linspace( 0, spectrum.shape[0], spectrum.shape[0] ) )
    
    return time, wavelenght, spectrum
    
    

def getM1600_old(_x0, spectremodeleenergparsparAngstrom):
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


def getModel(BPASS=False, oldSB99=False, filename=''):
    
    if(BPASS and not(oldSB99)):
        age, wavelength, spectre = readBPASS(filename=filename)
        
    if(not(BPASS) and oldSB99):
        age, wavelength, spectre = readStarburst_old(filename=filename)
        
    if(not(BPASS) and not(oldSB99)):
        age, wavelength, spectre = readStarburst(filename=filename)
      
    if(BPASS and oldSB99):
        print( "choose BPASS OR SB99 file!" )
        return -1

def getM1600(_x0, spectremodeleenergparsparAngstrom):
    lmin=1500
    lmax=1600
    return getMagInFilter(_x0, spectremodeleenergparsparAngstrom, lmin, lmax)



def getMagFilter(filter):
    """
    return the magnitude in a color band

    filter from SDSS
    source : http://www.cfht.hawaii.edu/Science/mswg/filters.html
    """

    if filter == "u":
        lmean = 3540
        dl=570./2
    if filter == "g":
        lmean = 4770
        dl=1370./2
    if filter == "r":
        lmean = 6230
        dl=1370./2
    if filter == "i":
        lmean = 7630
        dl=1530./2
    if filter == "z":
        lmean = 9130
        dl=950./2

    lmin=lmean-dl
    lmax=lmean+dl

    return getMagInFilter(_x0, spectremodeleenergparsparAngstrom, lmin, lmax)


def getMagInFilter(_x0, spectremodeleenergparsparAngstrom, lmin, lmax):
    parsec=3.08567758e16#parsec in meter
    c=299792458 #light speed

    F0 = 3631 #Jy reference flow
    Jyhzcm=1.e-23 #  erg s-1 Hz-1 cm-2
    Jyhzm=1.e-19 # erg s-1 Hz-1 m-2
    #so for absolute mags the 0 point is Jyhzm * 4pir^2 where r = 10 pc

    AB0pointabs=F0*Jyhzm*4.*np.pi*(10.*parsec)**2; # erg/s/Hz // cause ParSec

    AB0pointref=AB0pointabs*(c*1.e10/((_x0)**2)); # in erg/s/A (thats why we have to put c in Angstrom/s)

    mask = np.where( (_x0>=lmin) & (_x0<=lmax))#le filtre vaut 1 entre 1500 et 1600 Angstrom

    mag=-2.5*np.log10(integrate.trapz(spectremodeleenergparsparAngstrom[mask],
                                          _x0[mask])/integrate.trapz(AB0pointref[mask],
                                          _x0[mask]))
    return mag

def getModel2():
	""" LOOK HERE """
    age, wavelength, spectre = readStarburst()

    n = len(age)
    modelMag = np.zeros(n)

    for i in range(n):
        modelMag[i] = getM1600(wavelength, spectre[i])

    return  modelMag, age

def sumMag(mass, age, modelmag, modelage):
    mags=np.interp(age,modelage,modelmag)
    res=np.sum(mass*np.power(10,-mags/2.5))
    return -2.5*np.log10(res)

def get_all_flux_1600(stars,current_time,unit_mass, model='' ):
    """
    get luminous flux of stars in M1600 band
    """
    if( model=='' ):
        modelmag, modelage = getModel()
    else:
        modelmag, modelage = model
    age = current_time - stars.age.data
    mass = stars.mass.data*unit_mass/1.9891e30/1e6
    mags = np.interp(age,modelage,modelmag)
    flux = mass * mag2flux(mags)

    stars.flux_1600=flux

def get_all_flux_UV(mass,age,current_time,unit_mass,run):
    """
    get luminous flux of stars in UV
    """

    def getflux_UV(age,tlife,E0):
        y=np.ones(len(age)) *E0
        y[age>tlife] *= np.power(age[age>tlife]/tlife ,-4.)
        return y

    tlife = run.tlife_rad #year
    E0 = run.src_int_or_fesc*run.fesc #phot/s/kg

    mass *= unit_mass # kg
    flux = getflux_UV(current_time-age,tlife,E0) # #phot/s/kg
    return  flux*mass #phot/s

def get_all_flux_UV_old(stars,current_time,unit_mass,run):
    """
    get luminous flux of stars in UV
    """

    def getflux_UV(age,tlife,E0):
        y=np.ones(len(age)) *E0
        y[age>tlife] *= np.power(age[age>tlife]/tlife ,-4.)
        return y

    tlife = run.tlife_rad #year
    E0 = run.src_int_or_fesc #phot/s/kg

    mass = stars.mass.data*unit_mass # kg
    flux = getflux_UV(current_time - stars.age.data,tlife,E0) # #phot/s/kg
    stars.flux_UV = mass * flux #phot/s

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
