# coding: utf-8

import numpy as np
from scipy import integrate

class Constantes :
    """
    some physical constants
    """
    def __init__(self):

        self.c      = 299792458.       # velocity of light in m/s
        self.planck = 6.62617e-34      # J/s Planck constant
        self.M0     = 1.9891e30        # Solar masse
        self.eV     = 1.602176565e-19  # electron Volt in Joules
        self.erg    = 1.e-7            # erg in Joules
        self.Ang    = 1.e-10           # Angstreom in meter

cst = Constantes()

def lambda2E(l):
    """
    convert wavelenght[m] to energy [J]
    """
    cst=Constantes()
    c=cst.c
    h=cst.planck
    return h*c/l

def eV2lambda(eV):
    """
    convert energy [eV] to wavelenght [m]
    """
    cst=Constantes()
    h = cst.planck
    c = cst.c
    Ej = eV*cst.eV
    return h*c/Ej

def cross_section(egy):
    """
    compute cross section [m2]
    """
    P=2.963
    x=egy/4.298e-1

    a=(x-1.)**2
    b=np.power(x,P/2.-5.5)
    c=np.power(1+np.sqrt(x/3.288e1),-P)

    return 5.475e4*a*b*c*1e-22*(egy >= 13.6) # m2

def mag2flux(mag):
    """
    Transform Magnitude to Flux
    """
    flux = np.power(10,-mag/2.5)
    return flux

def flux2mag(flux):
    """
    Transform Flux to Magnitude
    """
    mag = -2.5*np.log10(flux)
    return mag

class  SP_Model:
    """
    Star Population Model
    give a file in input
    """
    def __init__(self, filename='', BPASS=False, SB99_old=False ):
        
        self.filename = filename
        self.BPASS = BPASS ### read BPASS file
        self.SB99_old = SB99_old ### read "old" SB99 file = fig7e and fig1e
        
        ### read file
        ### init time, wavelenght, spectrum
        self._getModel(  )
            
        ### automatic init Mag_1600 (basic Mag for everybody :)
        self.METH_getModelMag(  )
        
        #self.METH_compute_ionizing_properties() ### EMMA RAD INPUT
    
    #######################
    ### PRIVATE methods ###
    #######################
    ### i.e. should not be used out of the class OR you are doing something nasty
    
    def _getModel( self ):
        """
        'Parser' for starPopulation files 
        (Est ce que 'Parser est le bon mot ?')
        """

        if( self.BPASS and not(self.SB99_old) ):
            age, wavelength, spectre = self._readBPASS(  )

        if(not(self.BPASS) and self.SB99_old):
            age, wavelength, spectre = self._readStarburst_old(  )

        if(not(self.BPASS) and not(self.SB99_old)):
            age, wavelength, spectre = self._readStarburst(  )

        if(self.BPASS and self.SB99_old):
            print("ERROR : CHOOSE BPASS or SB99_old")
            raise AttributeError
            return -1

        self.time = age              ### [yr]
        self.wavelenght = wavelength ### [A]
        self.spectrum = spectre      ### [erg.s-1.A-1]
        

    def _readStarburst_old( self, convert=False ):
        """
        Reader for Starburst99 fig7e.dat file

        time [y]
        wavelenght [A]
        spectrum [erg.s-1.A-1]
        """

        if( self.filename=='' ):
            self.filename="/astro/home/nicolas.gillet/WORK/analyse_simu/fig7e.dat"

        print( 'Reading %s'%filename )
        #reading header
        with open(self.filename, 'r') as file:
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

        return age, data[0], data[1:]

    def _readStarburst( self ):
        """
        Reader for Starburst99 TH_AUBERT.spectrum1 file
        Files return by a 'personnal' Starburst simulation

        time [y]
        wavelenght [A]
        spectrum [erg.s-1.A-1]
        """

        # file = "../jupyter/Starburst99/Salpeter/SALP_AUBERT.spectrum1"
        # file = "../jupyter/Starburst99/Topheavy/TH_AUBERT.spectrum1"
        # file = "../jupyter/Starburst99/Kroupa/KROUPA_DOM.spectrum1"
        if( self.filename=='' ):
            self.filename = "/astro/home/nicolas.gillet/WORK/analyse_simu/TH_AUBERT.spectrum1"

        print( 'Reading %s'%self.filename )

        data = np.loadtxt( self.filename, unpack=True, skiprows=6 )
        time = np.unique(data[0,:])

        spectrum = []
        for cur_time in time:
                mask=data[0,:]==cur_time
                cur_data = data[:,mask]
                wavelenght = cur_data[1]
                spectrum.append(np.power(10,cur_data[2]))

        return time, wavelenght, np.array(spectrum)

    def _readBPASS( self ):
        """
        read BPASS spectrum

        time [y]
        wavelenght [A]
        spectrum [erg.s-1.A-1]
        """

        if self.filename=='':
            self.filename = "/astro/home/nicolas.gillet/WORK/analyse_simu/BPASSv2_imf135_300/OUTPUT_POP/spectra-bin.z001.dat"

        print( 'Reading %s'%self.filename )

        data = np.loadtxt(self.filename, unpack=True)

        spectrum = data[1:,:] * 3.8274e26 * 1.e7 ### convert from Luminosity_sun.A-1 to erg.s-1.A-1
        wavelenght = data[0]
        time = 10**( 6 +0.1*np.linspace( 0, spectrum.shape[0], spectrum.shape[0] ) )

        return time, wavelenght, spectrum
        
    def _getMagFilter( self, wavelenght, spectrum, Mag_filter ):
        """
        return the magnitude in a color band

        filter from SDSS
        source : http://www.cfht.hawaii.edu/Science/mswg/filters.html
        """
        if Mag_filter == "u":
            lmean = 3540
            dl=570./2
        if Mag_filter == "g":
            lmean = 4770
            dl=1370./2
        if Mag_filter == "r":
            lmean = 6230
            dl=1370./2
        if Mag_filter == "i":
            lmean = 7630
            dl=1530./2
        if Mag_filter == "z":
            lmean = 9130
            dl=950./2
        if Mag_filter == "1600":
            lmean = 1550
            dl = 100./2

        lmin=lmean-dl
        lmax=lmean+dl
        
        return self._getMagInFilter( wavelenght, spectrum, lmin, lmax )

    def _getMagInFilter( self, wavelenght, spectrum, lmin, lmax ):
        """
        return the magnitude between [lmin;lmax]
        """
        parsec = 3.08567758e16 ### parsec in meter
        c = 299792458 ### light speed

        F0 = 3631 ### Jy reference flow
        Jyhzcm=1.e-23 ###  erg s-1 Hz-1 cm-2
        Jyhzm=1.e-19 ### erg s-1 Hz-1 m-2
        ### so for absolute mags the 0 point is Jyhzm * 4pir^2 where r = 10 pc

        AB0pointabs=F0*Jyhzm*4.*np.pi*(10.*parsec)**2; ### erg/s/Hz // cause ParSec

        AB0pointref=AB0pointabs*(c*1.e10/((wavelenght)**2)); ### in erg/s/A (thats why we have to put c in Angstrom/s)

        mask = np.where( (wavelenght>=lmin) & (wavelenght<=lmax)) ### le filtre vaut 1 entre 1500 et 1600 Angstrom

        mag=-2.5*np.log10(integrate.trapz( spectrum[mask],
                                           wavelenght[mask]) / integrate.trapz( AB0pointref[mask],
                                                                                wavelenght[mask]))
        return mag
        
    ####################
    ### OPEN methods ###
    ####################

    def METH_getModelMag( self, Mag_filter="1600" ):
        """
        return the Mag from the spectrum
        
        filter from SDSS
        source : http://www.cfht.hawaii.edu/Science/mswg/filters.html
        """
        modelMag = np.zeros( len(self.time) )
        for i, spect in enumerate(self.spectrum):
            modelMag[i] = self._getMagFilter( self.wavelenght, spect, Mag_filter=Mag_filter )
        
        setattr( self, 'Mag_'+ Mag_filter, modelMag )
        
    def METH_get_all_flux_1600( self, stars, current_time, unit_mass ):
        """
        get luminous flux of all stars in M1600 band
        USE to compute Mag1600 in EMMA simulations
        """   
        age  = current_time - stars.age.data
        mass = stars.mass.data*unit_mass/1.9891e30/1e6 ### /!\ ATTENTION si le model n'a pas une masse de 1.e6!
        mags = np.interp( age, self.time, self.Mag_1600 )
        flux = mass * mag2flux( mags )
        stars.flux_1600 = flux
        
    def METH_compute_ionizing_properties( self ):
        """        
        return photon energy and cross section
        EMMA RAD INPUT
        
        time [yr]
        wavelenght [A]
        spectrum [erg.s-1.A-1]
        """

        #get constant
        cst = Constantes()

        #######################################################
        ### Compute ionizing properties at each stellar age ###
        #######################################################

        ### Select the ionizing photons
        ### 13.6 H ionization
        ### (lambda2E(10*cst.Ang)/cst.eV) = maximum limit 10A => BPASS spectrum seems to diverge ?
        mask = ((lambda2E(self.wavelenght*cst.Ang)/cst.eV) >= 13.6) \
        * ((lambda2E(self.wavelenght*cst.Ang)/cst.eV) <= (lambda2E(1*cst.Ang)/cst.eV))
        
        self.mask_lambda = mask

        x = self.wavelenght[mask] ### A
        y = self.spectrum[:,mask] ### erg.s-1.A-1

        ### total ionising energy
        E_ion = integrate.trapz( y, x ) ### erg.s-1

        ### total ionising photons
        Nphot_ion = integrate.trapz( y / (lambda2E( x *cst.Ang)/cst.erg), x ) ### Nph.s-1

        ### average energy of a ionizing photon
        Ephot = E_ion / Nphot_ion * cst.erg / cst.eV ### eV

        ### number of ionising photons per Msun
        nphot = Nphot_ion / (1e6*cst.M0) ### N.s-1.kg-1
        print( 'Warning : model mass assumed to be 1e6 Msun!' )

        ### photoionisation cross section
        s_lambda = cross_section( lambda2E( x *cst.Ang )/cst.eV ) ### m2

        ### number of photons per wavelenght bins
        N_lambda = y *cst.erg / lambda2E( x *cst.Ang )  ### Nph.s-1.A-1
        ### compute sigma, number weighted
        s_N = integrate.trapz( s_lambda*N_lambda, x ) / Nphot_ion ### m2

        ### compute sigma, energy weighted
        s_E = integrate.trapz( s_lambda*y, x ) / E_ion ### m2


        ####################################################################
        ### Compute ionizing properties average on the stellar evolution ###
        ####################################################################

        ### total number of ionizing photons per Msun
        Nphot = integrate.trapz( nphot, self.time*(3600*24*365.) ) ### Nph.kg-1
        ### average photon energy
        E_tot = integrate.trapz( Ephot*nphot, self.time*(3600*24*365.) ) / Nphot ### eV
        ### average sigma, number weighted
        s_N_tot = integrate.trapz( s_N*nphot, self.time*(3600*24*365.) ) / Nphot
        ### average sigma, energy weighted
        s_E_tot = integrate.trapz( s_E*nphot, self.time*(3600*24*365.) ) / Nphot
        
        self.Nphot = Nphot
        self.E_tot = E_tot
        self.s_N_tot = s_N_tot
        self.s_E_tot = s_E_tot

        print('')
        print( "FIRST age %.2f"%(self.time[0]/1.e6) )
        print( "hnu            Se                 Si")
        print("%s  %s  %s"%( Ephot[0], s_N[0], s_E[0]))
        print('')
        print( "average on all ages" )
        print( "hnu            Se                 Si")
        print("%s  %s  %s"%( E_tot, s_N_tot, s_E_tot))

        self.nphot = nphot
        self.Ephot = Ephot
        self.s_N = s_N
        self.s_E = s_E