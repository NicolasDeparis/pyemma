# coding: utf-8

import  numpy as np
import matplotlib.pyplot as plt

def sfr1():
    """
    data from Bouwens et al. 2015 http://arxiv.org/pdf/1403.4295v4.pdf

    format:
        z, log10(sfr)_dust_uncorrected, error_inf, error_sup, log10(sfr)_dust_corrected, error_inf, error_sup
    """
    z=np.array([
    [3.8,	-1.38, 	0.06, 0.06,	-1.00, 	0.06,	0.06],
    [4.9,  	-1.60,	0.06, 0.06,	-1.26, 	0.06,	0.06],
    [5.9, 	-1.80,	0.06, 0.06,	-1.55, 	0.06,	0.06],
    [6.8,  	-1.92,	0.06, 0.06,	-1.69, 	0.06,	0.06],
    [7.9,  	-2.23,	0.07, 0.06,	-2.08, 	0.07,	0.07],
    [10.4, 	-3.28,	0.36, 0.45, -3.13,	0.36,	0.45],
    # [10.5, 	-2.5,	0.5, 1. , -2.5,	0.5,	0.5]
    ])


    x=z[:,0]

    y1=np.power(10., z[:,1])
    low_err1 = np.abs(y1 - np.power(10., z[:,1] - z[:,3]))
    sup_err1 = np.abs(y1 - np.power(10., z[:,1] + z[:,2]))
    yerror1 = [low_err1, sup_err1]

    y2=np.power(10., z[:,4])
    low_err2 = np.abs(y2 - np.power(10., z[:,1] - z[:,6]))
    sup_err2 = np.abs(y2 - np.power(10., z[:,1] + z[:,5]))
    yerror2 = [low_err2, sup_err2]

    #plt.errorbar(x, y, yerr=yerror, ls='none', fmt='b*',label = "Bouwens et al. 2015")
    #plt.errorbar(x,y, yerr=yerror, ls='none', fmt='b*')

    x = np.append(x,x)
    y = np.append(y1,y2)

    low_err = np.append(low_err1,low_err2)
    sup_err = np.append(sup_err1,sup_err2)
    yerror = [low_err, sup_err]


#     print (x,y, yerror)

    return x, y, None, yerror

def sfr2():
    """
    Observational constraints for SFR.
    data from  Madau, P., & Dickinson, M. 2014, ARA&A, 52, 415 (MD 14)

    format:
        zmin, zmax, log10(sfr), error_inf, error_sup
    """
    z=np.array([
    [0.01,	0.1, 	-1.82,	+0.09,	-0.02],
    [0.2,	0.4,  	-1.50,	+0.05,	-0.05],
    [0.4,	0.6,  	-1.39,	+0.15,	-0.08],
    [0.6,	0.8,  	-1.20,	+0.31,	-0.13],
    [0.8,	1.2,  	-1.25,	+0.31,	-0.13],
    [0.05,	0.05,	-1.77,	+0.08,	-0.09],
    [0.05,	0.2,  	-1.75,	+0.18,	-0.18],
    [0.2,	0.4,  	-1.55,	+0.12,	-0.12],
    [0.4,	0.6,  	-1.44,	+0.10,	-0.10],
    [0.6,	0.8,  	-1.24,	+0.10,	-0.10],
    [0.8,	1.0,  	-0.99,	+0.09,	-0.08],
    [1.0,	1.2,  	-0.94,	+0.09,	-0.09],
    [1.2,	1.7,  	-0.95,	+0.15,	-0.08],
    [1.7,	2.5,  	-0.75,	+0.49,	-0.09],
    [2.5,	3.5,  	-1.04,	+0.26,	-0.15],
    [3.5,	4.5,  	-1.69,	+0.22,	-0.32],
    [0.92,	1.33, 	-1.02,	+0.08,	-0.08],
    [1.62,	1.88, 	-0.75,	+0.12,	-0.12],
    [2.08,	2.37, 	-0.87,	+0.09,	-0.09],
    [1.9,	2.7, 	-0.75,	+0.09,	-0.11],
    [2.7,	3.4, 	-0.97,	+0.11,	-0.15],
    [3.8,	3.8, 	-1.29,	+0.05,	-0.05],
    [4.9,	4.9, 	-1.42,	+0.06,	-0.06],
    [5.9,	5.9, 	-1.65,	+0.08,	-0.08],
    [7.0,	7.0, 	-1.79,	+0.10,	-0.10],
    [7.9,	7.9, 	-2.09,	+0.11,	-0.11],
    [7.0,	7.0, 	-2.00,	+0.10,	-0.11],
    [8.0,	8.0, 	-2.21,	+0.14,	-0.14],
    [0.03, 	0.03, 	-1.72,	+0.02,	-0.03],
    [0.03,	0.03, 	-1.95,	+0.20,	-0.20],
    [0.40,	0.70,  	-1.34,	+0.22,	-0.11],
    [0.70,	1.00,  	-0.96,	+0.15,	-0.19],
    [1.00,	1.30,  	-0.89,	+0.27,	-0.21],
    [1.30,	1.80,  	-0.91,	+0.17,	-0.21],
    [1.80,	2.30,  	-0.89,	+0.21,	-0.25],
    [0.40,	0.70,  	-1.22,	+0.08,	-0.11],
    [0.70,	1.00,  	-1.10,	+0.10,	-0.13],
    [1.00,	1.30,  	-0.96,	+0.13,	-0.20],
    [1.30,	1.80,  	-0.94,	+0.13,	-0.18],
    [1.80,	2.30,  	-0.80,	+0.18,	-0.15],
    [0.00,	0.30,  	-1.64,	+0.09,	-0.11],
    [0.30,	0.45,  	-1.42,	+0.03,	-0.04],
    [0.45,	0.60,  	-1.32,	+0.05,	-0.05],
    [0.60,	0.80,  	-1.14,	+0.06,	-0.06],
    [0.80,	1.00,  	-0.94,	+0.05,	-0.06],
    [1.00,	1.20,  	-0.81,	+0.04,	-0.05],
    [1.20,	1.70,  	-0.84,	+0.04,	-0.04],
    [1.70,	2.00,  	-0.86,	+0.02,	-0.03],
    [2.00,	2.50,  	-0.91,	+0.09,	-0.12],
    [2.50,	3.00,  	-0.86,	+0.15,	-0.23],
    [3.00,	4.20,  	-1.36,	+0.23,	-0.50],
    ])


    x= (z[:,1]+z[:,0])/2.
    xerror = (z[:,1]-z[:,0])/2.

    y= np.power(10., z[:,2])

    low_err = np.abs(y - np.power(10., z[:,2] + z[:,4]))
    sup_err = np.abs(y - np.power(10., z[:,2] + z[:,3]))
    yerror = [low_err, sup_err]

    return x, y, xerror, yerror

    #plt.errorbar(x, y, xerr=xerror, yerr=yerror, ls='none', fmt='b8', label = "Madau et al. 2014")

def sfr():
    x1, y1, xerror1, yerror1 = sfr1()
    x2, y2, xerror2, yerror2 = sfr2()

    separe=0
    if separe:
        plt.errorbar(x1, y1, xerr=xerror1, yerr=yerror1, ls='none', fmt='b8', label = "Madau et al. 2014")
        plt.errorbar(x2, y2, xerr=xerror2, yerr=yerror2, ls='none', fmt='b8', label = "Bouwens et al. 2015")
    else:
        plt.errorbar(x1, y1, xerr=xerror1, yerr=yerror1, ls='none', fmt='k.', label = "Observations")
        plt.errorbar(x2, y2, xerr=xerror2, yerr=yerror2, ls='none', fmt='k.')

def roberson_2015_fit():
    """ from : Cosmic Reionization and Early Star-forming Galaxies: A Joint Analysis of New Constraints from Planck and the Hubble Space Telescope"""
    z=np.linspace(0,20,1000)

    """Plank case"""
    a=0.01376
    b=3.26
    c=2.59
    d=5.68
    y_fit= a*np.power(1+z, b)/(1+np.power((1+z)/c,d))
    plt.plot(z,y_fit)

    """WMAP case"""
    a=0.01306
    b=3.66
    c=2.28
    d=5.29
    y_fit= a*np.power(1+z, b)/(1+np.power((1+z)/c,d))
    #plt.plot(z,y_fit)

def roberson_2010_fit():
    """from : Early star-forming galaxies and the reionization of the Universe"""
    z=np.linspace(0,20,100)

    a = 0.009
    b = 0.27
    h = 2.5

    """metal poor case"""
    c = 3.7
    d = 7.4
    g = 1e-3
    y_fit  = (a + b*(z/c)**h ) /(1. + (z/c)**d ) + g
    plt.plot(z,y_fit)

    """metal rich case"""
    c = 3.4
    d = 8.3
    g = 1e-4
    y_fit  = (a + b*(z/c)**h ) /(1. + (z/c)**d ) + g
    plt.plot(z,y_fit)

def xion():
    """
    Observational constraints for xion.
    data from Bouwens 2015
    http://cdsads.u-strasbg.fr/abs/2015ApJ...811..140B
    """


    #er=[[0.5],[1]]
    #plt.errorbar(6., 0.01, xerr=er, ls='none', fmt='ko')
    #plt.errorbar(6., 0.5,  xerr=er, ls='none', fmt='ko')

    z=np.array([

    [5.03, 0.9999451, 0.0000142, -0.0000165],
    [5.25, 0.9999330, 0.0000207, -0.0000244],
    [5.45, 0.9999333, 0.0000247, -0.0000301],
    [5.65, 0.9999140, 0.0000365, -0.0000460],
    [5.85, 0.9998800, 0.0000408, -0.0000490],
    [6.10, 0.99957  , 0.00030  , -0.00030  ],
    [7.0 , 0.66     , 0.12     , -0.09 	   ]
    ])

    x=z[:,0]
    y=z[:,1]

    yerror = [-z[:,3], z[:,2]]

    return x, 1.-y, None, yerror

    # plt.errorbar(x, 1.-y, xerr=None, yerr=yerror, ls='none', fmt='ko', label="Fan 2006")


    """

    ### pour les limites
    taille_erreur = 0.2 ### taille de la fleche, en dex

    # 5.9 >0.89 Dark Gaps in Quasar Spectra McGreer et al. (2015)
    # 5.6 >0.91 Dark Gaps in Quasar Spectra McGreer et al. (2015)
    z2 = np.array([ [ 5.9, 0.89, 0.05, 0 ], [ 5.6, 0.91, 0.05, 0 ] ]) ### limites sup
    x2 = z2[:,0]
    y2 = z2[:,1]
    yerror2 = [ (1-z2[:,1])*(1-10**(-taille_erreur)), z2[:,3] ] ### (1-z2[:,1])*(1-10**(-taille_erreur)) => pour que toute les fleche est la même taille
    plt.errorbar(x2, 1.-y2, xerr=None, yerr=yerror2, ls='none', fmt='ko', uplims=True )

    # 6.24-6.42 <0.9 (2) Ly Damping Wing of Quasars Schroeder et al. (2013)
    # Higher-Redshift Constraints

    # Prevalence of Lyα Emission in Galaxies S14
    # 4. 8.0 QHII (z = 8) < 0.35 Prevalence of Lyα Emission in Galaxies S14
    z3 = np.array([ [ 6.3, 0.9, 0, 0.1 ], [ 8.0, 0.35, 0, 0.1 ] ]) ### limites inf
    x3 = z3[:,0]
    y3 = z3[:,1]
    yerror3 = [ z3[:,2], (1-z3[:,1])*(10**(taille_erreur)-1) ] ### (1-z3[:,1])*(10**(taille_erreur)-1) => pour que toute les fleche est la même taille
    plt.errorbar(x3, 1.-y3, xerr=None, yerr=yerror3, ls='none', fmt='ko', lolims=True )

    """


    """
    Continuity with Ionizing Emissivity Estimates at z = 4.75
    5. log10 N˙
    ion(z = 4.75) = 1050.99±0.45 s
    −1 Mpc−3 BB13
    Other Constraints on the Ionization History of the Universe Not Explicitly Usedb

    6.3 ≥0.5 Lyα Damping Wing of a GRB Totani et al. (2006)	McQuinn et al. (2008)
    6.6 ≥0.6 Lyα Emitters Ouchi et al. (2010)
    6.6 ≥0.5 Galaxy Clustering McQuinn et al. (2007),	Ouchi et al. (2010)

    7.0 0.32-0.64 Lyα-Emitter LFs Ota et al. (2008)
    7.0 ∼0.5 Prevalence of Lyα Emission in Galaxies Caruana et al. (2014)
    7.0 0.1-0.4 Prevalence of Lyα Emission in Galaxies Ono et al. (2012)
    7.0 <0.49 Prevalence of Lyα Emission in Galaxies P14
    7.0 <0.5 Prevalence of Lyα Emission in Galaxies R13c
    7.0 <0.5 Clustering of Lyα Emitting Galaxies Sobacchi & Mesinger (2015)
    7.1 ≤0.9 Near-Zone Quasar Mortlock et al. (2011),
    Bolton et al. (2011)
    8.0 <0.70 Prevalence of Lyα Emission in Galaxies Tilvi et al. (2014)
    a This table is a compilation of the constraints presented in the original papers under References, but with valuable guidance
    """


"""
#if z < zp 1 − QHII(z) ∝ (1 + z)
#if z ≥ zp QHII(z) ∝ exp(−λ(1 + z)). (2)

zp=6.1
l=0.73
qzp = 0.99986



N=100
Q=np.zeros(N)
i=0

z = np.linspace(4,14,N)
for curz in z:
    if curz<zp:
        Q[i]=1+ (1+curz)**3
    else:
        Q[i]=np.exp(-l*(1+curz))*30

    i+=1
plt.plot(z,Q)
"""

def luminosity_function(redshift):
    """
    observation LF in UV for Bouwens et al 2014
    """

    if redshift == 4:
        mag=[-22.69,-22.19,-21.69,-21.19,-20.69,-20.19,-19.69,-19.19,-18.69,-18.19,-17.69,-16.94,-15.94]
        LF=[0.000003,0.000015,0.000134,0.000393,0.000678,0.001696,0.002475,0.002984,0.005352,0.006865,0.010473,0.024580,0.025080]
        error=[0.000004,0.000009,0.000023,0.000040,0.000063,0.000113,0.000185,0.000255,0.000446,0.001043,0.002229,0.003500,0.007860 ]

    if redshift == 5:
        mag=[-23.11,-22.61,-22.11,-21.61,-21.11,-20.61,-20.11,-19.61,-19.11,-18.36,-17.36,-16.36]
        LF=[0.000002,0.000006,0.000034,0.000101,0.000265,0.000676,0.001029,0.001329,0.002085,0.004460,0.008600,0.024400 ]
        error=[0.000002,0.000003,0.000008,0.000014,0.000025,0.000046,0.000067,0.000094,0.000171,0.000540,0.001760,0.007160 ]

    if redshift == 6:
        mag=[-22.52,-22.02,-21.52,-21.02,-20.52,-20.02,-19.52,-18.77,-17.77,-16.77]
        LF=[0.000002,0.000015,0.000053,0.000176,0.000320,0.000698,0.001246,0.001900,0.006680,0.013640]
        error=[0.000002,0.000006,0.000012,0.000025,0.000041,0.000083,0.000137,0.000320,0.001380,0.004200]

    if redshift == 7:
        mag=[-22.16,-21.66,-21.16,-20.66,-20.16,-19.66,-19.16,-18.66,-17.91,-16.91]
        LF=[0.000001,0.000033,0.000048,0.000193,0.000309,0.000654,0.000907,0.001717,0.005840,0.008500]
        error=[0.000002,0.000009,0.000015,0.000034,0.000061,0.000100,0.000177,0.000478,0.001460,0.002940]

    if redshift == 8:
        mag=[-21.87,-21.37,-20.87,-20.37,-19.87,-19.37,-18.62,-17.62]
        LF=[0.000005,0.000013,0.000058,0.000060,0.000331,0.000533,0.001060,0.002740]
        error=[0.000003,0.000005,0.000015,0.000026,0.000104,0.000226,0.000340,0.001040]

    if redshift == 10:
        mag=[-21.23,-20.23,-18.23]
        LF=[0.000001,0.000010,0.000266]
        error=[0.000001,0.000005,0.000171]

    return mag,LF,error

def luminosity_function_fit(redshift):
    """
    observation LF in UV for Bouwens et al 2014
    """

    MB=np.linspace(-25,-10,num=128)

    if redshift== 5:
        ms=-21.17
        ps=0.74e-3
        a=-1.76
        LFB=ps*np.log(10)/2.5*10.**(-0.4*(MB-ms)*(a+1))*np.exp(-10**(-0.4*(MB-ms)))

    if redshift== 6:
        ms=-20.94
        ps=0.5e-3
        a=-1.87
        LFB=ps*np.log(10)/2.5*10.**(-0.4*(MB-ms)*(a+1))*np.exp(-10**(-0.4*(MB-ms)))

    if redshift== 9:
        ms=-20.45
        ps=0.1e-3
        a=-2.3
        LFB=ps*np.log(10)/2.5*10.**(-0.4*(MB-ms)*(a+1))*np.exp(-10**(-0.4*(MB-ms)))
    return MB, LFB

def optical_depth(get=0):
    """
    plot data if get==0
    plot data and return values if get==1
    """

    zbecker2013=[2.15 ,2.25 ,2.35,2.45,2.55,2.65,2.75,2.85,2.95,3.05,3.15,3.25 ,3.35,3.45,3.55,3.65,3.75,3.85,3.95,4.05,4.15,4.25,4.35,4.45,4.55,4.65,4.75,4.85]
    meanFbecker2013=[0.8806, 0.8590, 0.8304, 0.7968, 0.7810, 0.7545, 0.7371, 0.7167, 0.6966, 0.6670, 0.6385, 0.6031, 0.5762, 0.5548, 0.5325, 0.4992, 0.4723, 0.4470, 0.4255, 0.4030, 0.3744, 0.3593,0.3441, 0.3216,0.3009, 0.2881, 0.2419, 0.2225]
    errmeanFbecker2013=[0.0103,0.0098,0.0093,0.0089,0.0090,0.0088,0.0088,0.0086,0.0084,0.0082,0.0080,0.0079,0.0074,0.0071,0.0071,0.0069,0.0068,0.0072,0.0071,0.0071,0.0074,0.0075,0.0102,0.0094,0.0104,0.0117,0.0201,0.0151]
    teffbecker2013=-np.log(meanFbecker2013)
    errtaueffbecker2013=-np.log(np.array(meanFbecker2013)+np.array(errmeanFbecker2013))-np.array(teffbecker2013)

    zabs_beck_2014=[5.902,5.737,5.577,5.423,5.948,5.781,5.620,5.464,5.635,5.479,5.328,5.183,5.586,5.432,
    5.283,5.139,5.577,5.423,5.274,5.130,5.478,5.327,5.182,5.423,5.275,5.131,4.992,4.858,5.138,4.999,4.864,4.734,4.608,
    5.043,4.907,4.776,4.648,4.525,4.996,4.861,4.731,4.605,4.484,4.929,4.797,4.669,4.545,4.425,4.824,4.696,4.571,4.450,4.333,
    4.720,4.594,4.473,4.355,4.242,4.691,4.567,4.446,4.329,4.216,4.634,4.511,4.393,4.278,4.166,4.625,4.502,4.384,4.269,4.158,
    4.453,4.336,4.223,4.113,4.007,4.358,4.244,4.134,4.027,3.923,4.349,4.235,4.125,4.018,3.915,4.282,4.170,4.062,3.958,3.856,
    4.282,4.170,4.062,3.958,3.856,4.253,4.143,4.035,3.932,3.831,4.244,4.134,4.027,3.923,3.822]

    Fmoy_beck_2014=[0.01149, 0.01578, 0.03836, 0.10284, 0.00820, 0.01404, 0.02544, 0.07504,0.03103, 0.02652,
    0.05167, 0.07156, 0.09962, 0.07989, 0.09959, 0.16812, 0.02299, 0.04055, 0.06021, 0.06568, 0.04605,
    0.06339, 0.13982, 0.04361, 0.06883, 0.11174, 0.07700, 0.14021, 0.23508, 0.16729, 0.17614, 0.19190, 0.19291, 0.12661,
    0.16273, 0.15605, 0.25978, 0.31753, 0.25503, 0.12560,0.22512, 0.30364, 0.29008, 0.17196, 0.10740, 0.18277, 0.16960,
    0.37799, 0.19396, 0.18355, 0.35492, 0.28457,0.32831,0.15083,0.28959,0.28650,0.35385,0.31313,0.29507,0.23270,0.42480,
    0.38612,0.34395,0.27845,0.29665,0.22136,0.29759,0.34961,0.20429,0.20725,0.45577,0.33738,0.29428,0.27003,0.32692,
    0.40195,0.39365,0.47038,0.39491,0.42531,0.35357,0.34943,0.38447,0.45173,0.28621,0.46655,0.27414,0.45430,0.38280,
    0.27687,0.32834,0.45227,0.51813,0.46323,0.38739,0.42680,0.46304,0.41029,0.47062,0.29638,0.41984,0.40853,0.48802,
    0.29407,0.41240,0.38504,0.45380,0.47881]

    errFmoy_beck_2014=[0.00263,0.00209,0.00187,0.00178,0.00075,0.00056,0.00062,0.00067,0.00197,0.00206,0.00150,
    0.00136,0.00253,0.00227,0.00190,0.00196,0.00281,0.00270,0.00210,0.00249,0.00043,0.00043,0.00040,0.00674,
    0.00505,0.00586,0.00758,0.00519,0.00062,0.00069,0.00063,0.00071,0.00074,0.00058,0.00059,0.00059,0.00068,0.00059,0.00107,
    0.00083,0.00100,0.00101,0.00105,0.00120,0.00103,0.00124,0.00132,0.00109,0.00147,0.00157,0.00172,0.00161,0.00153,0.00085,
    0.00088,0.00085,0.00095,0.00096,0.00045,0.00040,0.00037,0.00036,0.00036,0.00190,0.00181,0.00196,0.00206,0.00235,0.00184,
    0.00185,0.00221,0.00222,0.00254,0.00090,0.00205,0.00216,0.00197,0.00191,0.00075,0.00087,0.00090,0.00091,0.00097,0.00042,
    0.00043,0.00051,0.00045,0.00053,0.00171,0.00179,0.00184,0.00200,0.00205,0.00244,0.00287,0.00301,0.00323,0.00329,0.00093,
    0.00104,0.00111,0.00116,0.00127,0.00070,0.00075,0.00073,0.00075,0.00079]

    errtaueff_beck_2014=-np.log(np.array(Fmoy_beck_2014)+np.array(errFmoy_beck_2014))--np.log(np.array(Fmoy_beck_2014))

    # 4 errorbars in becker 2014
    zabs_beck_2014_2=[6.0737-0.01,5.797,5.7466 + 0.02,5.7369 - 0.01]
    errFmoy_min=-np.log(2*np.array([0.00196,0.00307,0.00280,0.00312]))
    errFmoy_max=-np.log(np.array([0.0005,0.0021,0.0010,0.0018]))
    y = np.array([1,1,1,1])
    ytop = errFmoy_max-y
    ybot = y-errFmoy_min

    # 3 arrows in becker 2014 ULAS J0148+0600 (2 first) and SDSS J2315−0023 (last entry)
    zulasj01480600=[5.796,5.634,5.965]
    errFmoy_j01480600=[0.00037,0.00051,0.00273]
    dxarrow=[0,0,0]
    dyarrow=[0.25,0.25,0.25]

    z_abs_fan2006=[5.58,5.43,5.28,5.13,4.98,5.64,5.49,5.34,5.19,5.04,5.81,5.66,5.51,5.36,5.21,5.06,5.52,5.37,5.22,5.07,4.92,5.66,
    5.51,5.36,5.21,5.06,5.61,5.46,5.31,5.16,5.01,6.10,5.95,5.80,5.65,5.50,5.35,5.55,5.40,5.25,5.10,4.95,5.68,5.53,5.38,5.83,5.68,5.53,5.38,5.23,
    5.08,6.25,6.10,5.95,5.80,5.65,5.50,5.90,5.75,5.60,5.45,5.30,5.77,5.62,5.47,5.32, 5.17,5.73,5.58,5.43,5.28,5.13,5.71,5.56,5.41,5.26,5.11,
    5.66,5.51,5.36,5.21,5.06,5.85,5.70,5.55,5.40,5.25,5.93,5.78,5.63,5.48,5.33,5.62,5.47,5.32,5.17]

    Fmoy_fan2006=[0.0170,0.0573 ,0.0205 ,0.1243 ,0.1002 ,0.0823 ,0.0718 ,0.0961 ,0.0578 ,0.1567 ,0.0216  ,0.0440 ,0.0984 ,0.1192 ,
    0.0884 ,0.1285 ,0.0907 ,0.0348 ,0.0606  ,0.0751 ,0.1276 ,0.0883  ,0.1127  ,0.1661 ,0.1191 ,0.1765 ,0.0884,0.1041 ,0.0596,0.1165,0.1268  ,0.0012  ,
    0.0060 ,0.0260 ,0.0462 ,0.0661 , 0.1147,0.0686 ,0.0520 , 0.0427,0.0898 ,0.1139 ,0.0117 ,0.0519 ,0.0736  ,0.0116 ,0.1010  ,0.0742 ,0.1341 ,0.1323 ,0.0530 ,
    0.0015  ,0.0051 ,0.0038  ,0.0186 ,0.0433 ,0.0278  ,0.0108  ,0.0055 ,0.0248 ,0.0077 ,0.0776  ,0.0645 ,0.0690  ,0.0991 ,0.0864 ,0.1156,0.0224,0.0445,
    0.1215,0.1217,0.1293,0.0322,0.0665,0.0858,0.0690,0.1650,0.0714,0.0775,0.0895,0.1292,0.1509,0.0687,0.0729,0.0795,0.0802,0.0934,0.0125,0.0071,
    0.0402,0.0407,0.0546,0.0495,0.1015,0.1376,0.0869]

    if get==0:
        plt.scatter(z_abs_fan2006,-np.log(Fmoy_fan2006), marker='o', facecolors='none', edgecolors='k', label=r"$\mathrm{Fan \, 2006}$")
        plt.errorbar(zbecker2013,teffbecker2013,yerr=errtaueffbecker2013, fmt='go',mfc='None', mec='g', label=r"$\mathrm{Becker \, 2013}$")
        plt.errorbar(zabs_beck_2014,-np.log(Fmoy_beck_2014),yerr=errtaueff_beck_2014, fmt='bs' ,mfc='None', mec='b', label=r"$\mathrm{Becker \, 2015}$")
        plt.errorbar(zabs_beck_2014_2, y, yerr=(ybot, ytop),fmt='bs',mfc='None', mec='w')

        for i in range(len(zulasj01480600)):
            plt.arrow(zulasj01480600[i], -np.log(2*errFmoy_j01480600[i]),dxarrow[i],dyarrow[i], fc="b", ec="b",head_width=0.04, head_length=0.2)

    if get==1:
        z_all = np.concatenate((z_abs_fan2006,zbecker2013,zabs_beck_2014))
        t_all = np.concatenate((-np.log(Fmoy_fan2006),teffbecker2013,-np.log(Fmoy_beck_2014)))
        return(z_all, t_all)


def baryonic_fraction(info):
    """
    from Okamoto 2008 : http://mnras.oxfordjournals.org/content/390/3/920.full.pdf
    """

    TODO

    # fb = info.ob/info.om
    #
    # alpha = 2
    #
    # res= (1+ (2**(alpha/3)-1) * (M/Mc)**(-alpha)  )**(-3/alpha)



def stellar_mass_function(redshift):
    """
    from song et al 2015
    http://arxiv.org/pdf/1507.05636v1.pdf
    available redshift are : 4, 5, 6, 7, 8
    """

    z = [4, 5, 6, 7, 8]
    M = [7.25,7.75,8.25,8.75,9.25,9.75,10.25,10.75,11.25]

    data = [
    [
    [-1.57, +0.21, -0.16],
    [-1.47, +0.24, -0.21],
    [-1.47, +0.35, -0.32],
    [-1.63, +0.54, -0.54],
    [-1.73, +1.01, -0.84],
    ],[
    [-1.77, +0.15, -0.14],
    [-1.72, +0.20, -0.20],
    [-1.81, +0.23, -0.28],
    [-2.07, +0.45, -0.41],
    [-2.28, +0.84, -0.64],
    ],[
    [-2.00, +0.13, -0.10],
    [-2.01, +0.16, -0.16],
    [-2.26, +0.21, -0.16],
    [-2.49, +0.38, -0.32],
    [-2.88, +0.75, -0.57],
    ],[
    [-2.22, +0.09, -0.09],
    [-2.33, +0.15, -0.10],
    [-2.65, +0.15, -0.15],
    [-2.96, +0.32, -0.30],
    [-3.45, +0.57, -0.60],
    ],[
    [-2.52, +0.09, -0.09],
    [-2.68, +0.07, -0.14],
    [-3.14, +0.12, -0.11],
    [-3.47, +0.32, -0.35],
    [-4.21, +0.63, -0.78],
    ],[
    [-2.91, +0.12, -0.05],
    [-3.12, +0.09, -0.11],
    [-3.69, +0.12, -0.13],
    [-4.11, +0.41, -0.57],
    [-5.31, +1.01, -1.64],
    ],[
    [-3.41, +0.13, -0.08],
    [-3.63, +0.13, -0.11],
    [-4.55, +0.19, -0.24],
    [-5.12, +0.68, -1.10],
    [-6.93, +1.61, -2.57],
    ],[
    [-4.11, +0.22, -0.21],
    [-4.40, +0.15, -0.35],
    [-5.96, +0.52, -0.32],
    [-6.57, +1.14, -1.37],
    [-8.83, +2.57, -5.75],
    ],[
    [-5.00, +0.24, -0.97],
    [-5.96, +0.98, -1.98],
    [-7.50, +1.30, -0.99],
    [-8.59, +2.34, -3.23],
    [-12.10, +4.26, -13.54]
    ]]


    # data=np.power(10,data)

    cur_z=z.index(redshift)

    n=len(M)
    M_out = np.zeros(n)
    E_inf = np.zeros(n)
    E_sup = np.zeros(n)

    for cur_M in range(n):

        M_out[cur_M]=data[int(cur_M)][int(cur_z)][0]
        E_inf[cur_M]=data[int(cur_M)][int(cur_z)][2]
        E_sup[cur_M]=data[int(cur_M)][int(cur_z)][1]

    return M,M_out,E_inf,E_sup

def stellar_mass_function_fit():
    """
    from song et al 2015
    http://arxiv.org/pdf/1507.05636v1.pdf

    available redshift are : 4, 5, 6, 7, 8
    """
    ### zcen(binMs) => à remplacer : les centres des bins de masse (stellaire)
    ### all song (mimi) et al (2015) fits,
    ### z, m_norm [Msun], alpha, phi_norm [Mpc-3]
    SMF_fits = np.array([[ 4, 10**(10.43), -1.53, 30.50e-5 ],
                         [ 5, 10**(10.46), -1.66, 13.80e-5 ],
                         [ 6, 10**(10.32), -1.94, 2.78e-5 ],
                         [ 7, 10**(10.42), -2.04, 0.65e-5 ],
                         [ 8, 10**(10.41), -2.40, 0.03e-5 ]])
    ### apply to the fits
    SMF_mods = np.array([ (SMF_fits[0,3]/SMF_fits[0,1])*(zcen(binMs)/SMF_fits[0,1])**SMF_fits[0,2]*np.exp(-(zcen(binMs)/SMF_fits[0,1]) )*(binMs[1:]-binMs[:-1]),
                          (SMF_fits[1,3]/SMF_fits[1,1])*(zcen(binMs)/SMF_fits[1,1])**SMF_fits[1,2]*np.exp(-(zcen(binMs)/SMF_fits[1,1]) )*(binMs[1:]-binMs[:-1]),
                          (SMF_fits[2,3]/SMF_fits[2,1])*(zcen(binMs)/SMF_fits[2,1])**SMF_fits[2,2]*np.exp(-(zcen(binMs)/SMF_fits[2,1]) )*(binMs[1:]-binMs[:-1]),
                          (SMF_fits[3,3]/SMF_fits[3,1])*(zcen(binMs)/SMF_fits[3,1])**SMF_fits[3,2]*np.exp(-(zcen(binMs)/SMF_fits[3,1]) )*(binMs[1:]-binMs[:-1]),
                          (SMF_fits[4,3]/SMF_fits[4,1])*(zcen(binMs)/SMF_fits[4,1])**SMF_fits[4,2]*np.exp(-(zcen(binMs)/SMF_fits[4,1]) )*(binMs[1:]-binMs[:-1]),
                        ])



def ionization_rate_gamma():
    """
    Observations of Gamma the photo-ionization rate.
    The observation are not well calibrate to each other (see Becker & Bolton 13 Figure.11)
    but the difference is small
    """
    ### Becker & Bolton 13
    obs_gamma = np.array( [ [ 2.40, 2.80, 3.20, 3.60, 4.00, 4.40, 4.75 ],
                           [ 0.015, -0.066, -0.103, -0.097, -0.072, -0.019, -0.029 ],
                           [ -0.146, -0.131, -0.121, -0.118, -0.117, -0.122, -0.147 ],
                           [ 0.132, 0.129, 0.130, 0.131, 0.135, 0.140, 0.156 ]] )
    ### Calverley 11
    obs_gamma2 = np.array( [ [ 5, 6 ],
                             [ -0.15, -0.84 ],
                             [ -0.16, -0.18 ],
                             [  0.16, 0.18 ]] )
    ### Wyithe and Bolton 11
    obs_gamma3 = np.array( [ [ 5, 6 ],
                             [ np.log10(0.47), np.log10(0.18) ],
                             [ -0.2, -0.09 ],
                             [ 0.3, 0.18 ]] )

    plt.plot( obs_gamma[0], 10**obs_gamma[1], 'ok')
    plt.errorbar( obs_gamma[0], 10**obs_gamma[1],
                 yerr=[10**obs_gamma[1] - 10**(obs_gamma[1]+obs_gamma[2]),
                       10**(obs_gamma[1]+obs_gamma[3]) - 10**obs_gamma[1] ],
                fmt='none', ecolor='k' )
    plt.plot( obs_gamma2[0]-0.02, 10**obs_gamma2[1], 'or')
    plt.errorbar( obs_gamma2[0]-0.02, 10**obs_gamma2[1],
                 yerr=[10**obs_gamma2[1] - 10**(obs_gamma2[1]+obs_gamma2[2]),
                       10**(obs_gamma2[1]+obs_gamma2[3]) - 10**obs_gamma2[1] ],
                fmt='none', ecolor='r' )
    plt.plot( obs_gamma3[0]+0.02, 10**obs_gamma3[1], 'ob')
    plt.errorbar( obs_gamma3[0]+0.02, 10**obs_gamma3[1],
                 yerr=[10**obs_gamma3[1] - 10**(obs_gamma3[1]+obs_gamma3[2]),
                       10**(obs_gamma3[1]+obs_gamma3[3]) - 10**obs_gamma3[1] ],
                fmt='none', ecolor='b' )
