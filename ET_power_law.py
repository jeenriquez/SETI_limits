import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy import units as u
from matplotlib.patches import Polygon
from scipy import stats


def calc_DishArea(d):
    """ Compute dish area
    d = dish diameter
    """
    return np.pi * (d/2)**2

def calc_BeamSize(d,v,verbose=False):
    """ Compute BeamSize
    d = dish diameter
    v = frequency
    """

    c = 2.998e8  #speed of light

    if verbose:
        print '\nBeam is: %f \n'%(1.22* (c/(d*v)) *57.2958)

    return (1.22* (c/(d*v)) *57.2958/2)**2*np.pi

def calc_SEFD(A, Tsys, eff=1.0):
    """ Calculate SEFD
    Tsys = system temperature
    A = collecting area
    Ae = effective collecting area
    eff = aperture efficency (0.0 to 1.0)
    """
    kb = 1.3806488e3  # 1.38064852e-23    Boltzmann constant
    Ae = A*eff
    return 2 * Tsys * kb / Ae

def calc_Sensitivity(m, nu, t, SEFD=0, Tsys=10,eff=1.0,A=100,npol=2.,narrow=True):
    """ Minimum detectable luminosity for narrowband emission
    Tsys = system temperature
    A = collecting area
    m = threshold, (e.g. 10)
    nu = channel bandwidth
    t = observing time
    narrow = True if signal is narrower than spectral resolution.
    """
    if not SEFD:
        sefd = calc_SEFD(A, Tsys, eff=eff)
    else:
        sefd = SEFD


    if narrow:
        sens = m * sefd * np.sqrt(nu/(npol*t))
    else:
        sens = m * sefd / np.sqrt(npol*nu*t)

    return sens

def calc_EIRP_min(d,Sens):
    """ Minimum detectable luminosity (EIRP) for narrowband emission.
    d = distance to target star []
    Sens = sensitivity of the obs (Jy)
    """

    #1 Jy = 1e-26 W/m2/Hz)

    return 4 * np.pi * d**2 * Sens *1e-26

def calc_gain(ll, d):
    """ Gain of a dish telescope
    ll = wavelength (lambda)
    d = dish diameter (m)
    """
    return (np.pi * d / ll)**2

def calc_NP_law(alpha,P):
    ''' Calcualtes the power law, given a the power exponend and an array of EIRP powers in W. Based on Gray&Mooley
    '''

#    No = 1e25**alpha   # Not sure how I got 1e25 for the transmitter power, but from the plot I get a lower number.
    No = 1.e21**alpha
    NP = No*(1./P**(alpha))
    NPn = NP*1e3/4e11   # Normalized to stars in MW (4e11), and BW (1e3 .. ).

    return NPn

def calc_NP_law2(alpha,P):
    ''' Calcualtes the power law, given a the power exponend and an array of EIRP powers in W. Based on BL-Lband.
    '''

    No = 706146012574.0**alpha / 304.48
    NP = No*(1./P**(alpha))

    return NP

def calc_NP_law3(P):
    ''' Calcualtes the power law, given a the power exponend and an array of EIRP powers in Watts. Based on fit from both above.
    '''

    E1 = 5.e+11
    S1  = 350

    E2 = 1.98792219e+21  #2.75879335e+21
    S2 = 7.14285714e+08

    # Solving for alpha
    alpha = np.log10(S2/S1) /np.log10(E2/E1)
    print 'The exponent (alpha) = ', alpha
    # Solving for No
    No = E1**alpha / S1

    NP = No*(1./P**(alpha))

    return NP


def ET_power_law(verbose=False):

    #---------------------------
    # Standardizing the sensitivity of telescopes by figuring out the max EIRP of a transmitter they could have found.
    # For this, I need the sensitivity of the observation, given SEFD and so on.

    #Standard distance (100 ly)
    dist_std = (100. * u.lyr.to('m'))
    #So decided not to use the distance above, it does makes sense if the original distance is shorter, but not if the distance is farther away (like in Siemion2013)

    #nomalization to L_AO = 2e13 W, fraq_freq = 1/2 and N_stars = 1k
    zeta_AO  =  1 #np.log10(1e3*.5)/ np.log10(2e13)
    zeta_AO  =  1e3*0.5/ 1e13

    #---------------------------
    #BL
    telescope = 'GBT'
    BL_stars= 692
    band = 660e6  #(1.1-1.2,1.34-1.9)
    freq_range_norm = (band/1.5e9)

    BL_SEFD = calc_SEFD(calc_DishArea(100), 20, eff=0.72) # 10 Jy (GBT)
    m=25.; nu =3.; t=300.
    Sens = calc_Sensitivity(m, nu,t,SEFD=BL_SEFD)

    dist_std = (50*3.26156 * u.lyr.to('m'))  # Max distance approx
    BL_EIRP = calc_EIRP_min(dist_std,Sens)
    BL_rarity = BL_stars*freq_range_norm

    iband = 800e6
    BL_speed = BL_SEFD**2*nu/iband

    BL_sky = BL_stars * calc_BeamSize(100,1.5e9)
    BL_DFM = BL_sky * band / Sens**(3/2.)

    if verbose:
        print 'SEFD (BL):', BL_SEFD
        print 'Sens (BL):', Sens
        print 'EIRP (BL):', BL_EIRP
        print 'BeamSize (BL):', calc_BeamSize(100,1.5e9)
        print 'Sky Coverage (BL):', BL_sky
        print 'CWTFM (BL):',  zeta_AO *(BL_EIRP) / (BL_rarity)
        print 'DFM (BL):', BL_DFM
        print '~o~'

    #----------
    # Gray & Mooley 2017

    telescope=['VLA']
    GM_stars= 1e12

    band =  np.array([1e6, 0.125e6])
    central_freq = np.array([1.4e9,8.4e9])
    freq_range_norm = (band/central_freq)

    SEFD = calc_SEFD(calc_DishArea(25),35, eff=0.45) / np.sqrt(27*26.)   #Perley 2009
    m=7.0; nu =[122., 15.3]; t=[20*60,5*60]
    Sens = np.array([calc_Sensitivity(m, nu[0],t[0],SEFD=SEFD),calc_Sensitivity(m, nu[1],t[1],SEFD=SEFD)])

    dist_std = (2.5e6 * u.lyr.to('m'))  # Max distance approx
    GM_EIRP = np.array([calc_EIRP_min(dist_std,Sen) for Sen in Sens])

    GM_rarity = GM_stars*freq_range_norm

    GM_rarity_tot = GM_rarity.sum()
    GM_EIRP_tot = GM_EIRP.max()

    iband = 1e6
    GM_speed = SEFD**2*nu[0]/iband

    GM_sky = 8*(0.95/2.)**2*np.pi   # 0.95deg images    #NOTE: in Enriquez 2017, we used this as the radius, but it is the diameter.
    GM_DFM = GM_sky * band / Sens**(3/2.)

    if verbose:
        print 'SEFD (Gray & Mooley 2017):', SEFD
        print 'Sens (Gray & Mooley 2017):', Sens
        print 'EIRP (Gray & Mooley 2017):', GM_EIRP
        print 'BeamSize (Gray & Mooley 2017):',
        print 'Sky Coverage (Gray & Mooley 2017):', GM_sky
        print 'CWTFM (Gray & Mooley 2017):',  zeta_AO * (GM_EIRP_tot)/ (GM_rarity_tot) #,'or', zeta_AO*stats.hmean(GM_EIRP/GM_rarity)
        print 'DFM (Gray & Mooley 2017):', GM_DFM
        print '~o~'

    #----------
    #Phoenix
    telescope = ['Arecibo','Arecibo; Parkes,Parkes,NRAO140']

    Ph_stars = np.array([290,371,206,105,195]) # From Harp2016

    #180MHz skip Backus2002; band from Harp2016
    band = np.array([(1.75-1.2)*1e9 - 180e6,(3.0-1.75)*1e9,(1.75-1.2)*1e9, (3.0-1.75)*1e9, (3.0-1.2)*1e9])

    central_freq = np.array([1.5e9,2.375e9,1.5e9,2.375e9,2.1e9])
    freq_range_norm = (band/central_freq)

    Dish_D = np.array([305,225,64,64,43])  # Email from G. Harp

    SEFD = np.array([calc_SEFD(calc_DishArea(Dish_D[0]), 40, eff=0.7),
        calc_SEFD(calc_DishArea(Dish_D[1]), 40, eff=0.7),
        calc_SEFD(calc_DishArea(Dish_D[2]), 35, eff=0.7),
        calc_SEFD(calc_DishArea(Dish_D[3]), 35, eff=0.7),
        calc_SEFD(calc_DishArea(Dish_D[4]), 35, eff=0.7)])

    m=1; nu =1.0; t=[276,195,276,138,552]
    Sens1 = np.array([calc_Sensitivity(m,nu,t[i],SEFD=SEFD[i],narrow=False) for i in range(len(SEFD))])
    Sens = np.array([16,16,100,100,100])  # From Harp2016

    # Max distance approx    ; 147Ly median distance Shostalk(2000), ~700 farthest
    dist_std = (700 * u.lyr.to('m'))

    Ph_EIRP = np.array([calc_EIRP_min(dist_std,Sen) for Sen in Sens])
    Ph_rarity = Ph_stars*freq_range_norm
    Ph_stars_tot = Ph_stars.sum()
    Ph_rarity_tot = Ph_rarity.sum()
    Ph_EIRP_tot = Ph_EIRP.max()

    iband = 20e6
    Ph_speed = SEFD.mean()**2*nu/iband   #Note: This value is calculated with self calculated SEFD values (which are not completely consistent with values expected from Harp 2016 values).

    Ph_sky = Ph_stars * np.array([calc_BeamSize(Dish_D[i],central_freq[i]) for i in range(len(Dish_D))])
    Ph_DFM = Ph_sky * band / Sens**(3/2.)

    if verbose:
        print 'SEFD (Phoenix):', SEFD
        print 'Sens (Phoenix):', Sens1
        print 'Sens_Harp (Phoenix):', Sens
        print 'EIRP (Phoenix):', Ph_EIRP
        print 'BeamSize (Phoenix):', np.array([calc_BeamSize(Dish_D[i],central_freq[i]) for i in range(len(Dish_D))])
        print 'Sky Coverage (Phoenix):', Ph_sky.sum()
        print 'CWTFM (Phoenix):',  zeta_AO * (Ph_EIRP)/ (Ph_rarity)

        print 'CWTFM (Phoenix):',  zeta_AO * (Ph_EIRP_tot)/ (Ph_rarity_tot)
        print 'DFM (Phoenix):', Ph_DFM
        print '~o~'

    #----------
    #ATA
    telescope = 'ATA'

    ATA_stars= np.array([65,1959,2822,7459])

    band = np.array([8000.e6,2040.e6,337.e6,268.e6])    #There are 73MHz which are RFI flagged, it is ignored here.
    central_freq = 5e9    # 1-9 GHz
    freq_range_norm = (band/central_freq)

    #Tsys = (80+120+95+137)/4. = 108
    SEFD = calc_SEFD(calc_DishArea(6.1), 108, eff=0.58) / np.sqrt(27*26)
    SEFDs = np.array([SEFD,SEFD,SEFD,SEFD])

    m=6.5; nu =0.7; t=93.
    dist_std = np.array([(1.4e3*3.26156 * u.lyr.to('m')),(1.1e3*3.26156 * u.lyr.to('m')),(300 * u.lyr.to('m')),(500 * u.lyr.to('m'))])  #Turnbull 2003 for HabCat
    Sens = np.array([calc_Sensitivity(m,nu,t,SEFD=SEF,narrow=False) for SEF in SEFDs])
    ATA_EIRP = np.array([calc_EIRP_min(dist_std[i],Sens[i]) for i in range(len(Sens))])

    ATA_rarity = ATA_stars*freq_range_norm
    ATA_rarity_tot = ATA_rarity.sum()
    ATA_stars_tot = ATA_stars.sum()
    ATA_EIRP_tot = ATA_EIRP.max()

    iband = 70e6
    ATA_speed = SEFD**2*nu/iband

    ATA_sky = ATA_stars * 3*6./4.*np.pi/3600.  # beam 3'x6' at 1.4GHz
    ATA_DFM = ATA_sky * band / Sens**(3/2.)

    if verbose:
        print 'SEFD (ATA):', SEFD
        print 'Sens (ATA):', Sens
        print 'EIRP (ATA):', ATA_EIRP
        print 'BeamSize (ATA):',
        print 'Sky Coverage (ATA):', ATA_sky.sum()
        print 'CWTFM (ATA):',  zeta_AO * (ATA_EIRP_tot)/ (ATA_rarity_tot)
        print 'DFM (ATA):', ATA_DFM
        print '~o~'

    #----------
    #Siemion 2013
    telescope = 'GBT'
    Siemion_stars= 86

    band = 800e6 - 130e6  #(1.1-1.2,1.33-1.9)
    freq_range_norm = (band/1.5e9)

    SEFD= calc_SEFD(calc_DishArea(100), 20, eff=0.72) # 10 Jy (GBT)
    m=25.; nu =1.; t=300.
    Sens = calc_Sensitivity(m, nu,t,SEFD=SEFD)

    dist_std = (1.1e3*3.26156 * u.lyr.to('m'))  # Max distance approx
    Siemion_EIRP = calc_EIRP_min(dist_std,Sens)

    Siemion_rarity = Siemion_stars*freq_range_norm

    iband = 800e6
    Siemion_speed = (SEFD/0.85)**2*nu/iband    # 0.85 ==> Siemion priv. comm. (due to 2 bit data format)

    Siemion_sky = Siemion_stars * calc_BeamSize(100,1.5e9)
    Siemion_DFM = Siemion_sky * band / Sens**(3/2.)

    if verbose:
        print 'SEFD (Siemion2013):', SEFD
        print 'Sens (Siemion2013):', Sens
        print 'EIRP (Siemion2013):', Siemion_EIRP
        print 'BeamSize (Siemion2013):',calc_BeamSize(100,1.5e9)
        print 'Sky Coverage (Siemion2013):', Siemion_sky
        print 'CWTFM (Siemion2013):',  zeta_AO * (Siemion_EIRP)/ (Siemion_rarity)
        print 'DFM (Siemion2013):', Siemion_DFM
        print '~o~'

    #----------
    #Valdes 1986
    telescope='HCRO'
    Valdes_stars = np.array([53, 12])

    band = np.array([256*4883, 1024*76])
    freq_range_norm = (band/1.516e9)

    SEFD = calc_SEFD(calc_DishArea(26), 100, eff=0.5)
    m=3.0; nu =[4883., 76.]; t=3000.
    Sens = np.array([calc_Sensitivity(m, nu[0],t,SEFD=SEFD,npol=1.),calc_Sensitivity(m, nu[1],t,SEFD=SEFD,npol=1.)])

    dist_std = (20 * u.lyr.to('m'))  # Max distance approx
    Valdes_EIRP = np.array([calc_EIRP_min(dist_std,Sen) for Sen in Sens])
    Valdes_rarity = Valdes_stars*freq_range_norm

    Valdes_rarity_tot = Valdes_rarity.sum()
    Valdes_EIRP_tot = Valdes_EIRP.max()

    iband = 256*4883
    Valdes_speed = SEFD**2*nu[0]/iband

    Valdes_sky = (Valdes_stars * calc_BeamSize(26,1.5e9)).sum()
    Valdes_DFM = Valdes_sky * band / Sens**(3/2.)

    if verbose:
        print 'SEFD (Valdes 1986):', SEFD
        print 'Sens (Valdes 1986):', Sens
        print 'EIRP (Valdes 1986):', Valdes_EIRP
        print 'BeamSize (Valdes 1986):',calc_BeamSize(26,1.5e9)
        print 'Sky Coverage (Valdes 1986):', Valdes_sky
        print 'CWTFM (Valdes 1986):',  zeta_AO * (Valdes_EIRP_tot)/ (Valdes_rarity_tot)
        print 'DFM (Valdes 1986):', Valdes_DFM

        print '~o~'

    #----------
    #Tarter 1980
    telsecope = 'NRAO 91m'
    Tarter_stars=201

    band = 360e3*4.
    freq_range_norm = (band/1.666e9)

    SEFD = calc_SEFD(calc_DishArea(91), 70, eff=0.6)
    m=12.0; nu =5.5; t= 45
    Sens = calc_Sensitivity(m, nu,t,SEFD=SEFD)

    dist_std = (25*3.26156* u.lyr.to('m'))  # Max distance approx
    Tarter_EIRP = calc_EIRP_min(dist_std,Sens)
    Tarter_rarity = Tarter_stars*freq_range_norm

    iband = 360e3
    Tarter_speed = SEFD**2*nu/iband

    Tarter_sky = Tarter_stars * calc_BeamSize(91,1.666e9)
    Tarter_DFM = Tarter_sky * band / Sens**(3/2.)

    if verbose:
        print 'SEFD (Tarter1980):', SEFD
        print 'Sens (Tarter1980):', Sens
        print 'EIRP (Tarter1980):', Tarter_EIRP
        print 'BeamSize (Tarter1980):', calc_BeamSize(91,1.666e9)
        print 'Sky Coverage (Tarter1980):', Tarter_sky
        print 'CWTFM (Tarter1980):',  zeta_AO * (Tarter_EIRP)/ (Tarter_rarity)
        print 'DFM (Tarter1980):', Tarter_DFM

        print '~o~'

    #----------
    #Verschuur1973
    telescope=['300ft Telescope', '140ft Telescope']
    Verschuur_stars=np.array([3,8])

    band = np.array([0.6e6,20e6])
    freq_range_norm = (band/1.426e9)

    SEFD = np.array([calc_SEFD(calc_DishArea(300*0.3048),110, eff=0.75),calc_SEFD(calc_DishArea(140*0.3048),48, eff=0.75)]) #**NOTE** the 0.75 for the 140' is not real
    m=3.0; nu =[490.,7.2e3]; t= [4*60.,5*60.]
    Sens = np.array([calc_Sensitivity(m, nu[0],t[0],SEFD=SEFD[0]),calc_Sensitivity(m, nu[1],t[1],SEFD=SEFD[1])])

    dist_std = (5*3.26156 * u.lyr.to('m'))
    Verschuur_EIRP = np.array([calc_EIRP_min(dist_std,Sen) for Sen in Sens])

    Verschuur_rarity = Verschuur_stars*freq_range_norm

    Verschuur_rarity_tot = Verschuur_rarity.sum()
    Verschuur_EIRP_tot = Verschuur_EIRP.max()

    iband = np.array([0.6e6, 2.5e6])  #300 ft: Two 192-channel receivers (at 130 km/s with 4.74kHz=1km/s at this freq.)
    Verschuur_speed = SEFD.min()**2*nu[0]/iband[0]

    Verschuur_sky = (Verschuur_stars * np.array([calc_BeamSize(300*0.3048,1.42e9),calc_BeamSize(140*0.3048,1.42e9)])).sum()*2  # The two comes from the off beam.
    Verschuur_DFM = Verschuur_sky * band / Sens**(3/2.)

    if verbose:
        print 'SEFD (Verschuur1973):', SEFD
        print 'Sens (Verschuur1973):', Sens
        print 'EIRP (Verschuur1973):', Verschuur_EIRP
        print 'BeamSize (Verschuur1973):', np.array([calc_BeamSize(300*0.3048,1.42e9),calc_BeamSize(140*0.3048,1.42e9)])
        print 'Sky Coverage (Verschuur1973):', Verschuur_sky
        print 'CWTFM (Verschuur1973):',  zeta_AO * (Verschuur_EIRP_tot)/ (Verschuur_rarity_tot)
        print 'DFM (Verschuur1973):', Verschuur_DFM

        print '~o~'

    #----------
    #META Horowitz&Sagan

    telescope=''
    Horowitz_stars= 1e7

    band = 1.2e6
    freq_range_norm = (band/1.42e9)

    SEFD = calc_SEFD(calc_DishArea(26),85, eff=0.5)    # eff=0.5 ==>  We were unable to find a value in the literature. We assume a similar value to the antenna of the same dimensions from Valdes & Freitas (1986).
    m=30; nu =0.05; t=20
    Sens = calc_Sensitivity(m, nu,t,SEFD=SEFD,narrow=False)

    dist_std = (700*3.26156 * u.lyr.to('m'))  # Max distance: # Horowitz & Sagan (1993) suggested values for the number of stars given a distance, based on the power of an isotropic beacon.
    Horowitz_EIRP = calc_EIRP_min(dist_std,Sens)
    Horowitz_rarity = Horowitz_stars*freq_range_norm

    iband = 400e3
    Horowitz_speed = SEFD**2*nu/iband

    Horowitz_sky = 41253*.68 # Horowitz_stars * calc_BeamSize(26,1.42e9)
    Horowitz_DFM = Horowitz_sky * band / Sens**(3/2.)

    if verbose:
        print 'SEFD (Horowitz):', SEFD
        print 'Sens (Horowitz):', Sens
        print 'EIRP (Horowitz):', Horowitz_EIRP
        print 'BeamSize (Horowitz):',
        print 'Sky Coverage (Horowitz):', Horowitz_sky
        print 'CWTFM (Horowitz):',  zeta_AO * (Horowitz_EIRP)/ (Horowitz_rarity)
        print 'DFM (Horowitz):', Horowitz_DFM
        print '~o~'

    #---------------------------
    #BL
    Price_telescopes = ['GBT','GBT','Parkes']

    Price_BL_stars = np.array([883,1006,195])
    Price_band = np.array([(1.2-1.025+1.925-1.34)*1e9,(2.72-1.82)*1e9,(3.444-2.574)*1e9])

    Price_central_freq = np.array([1.5e9,2.27e9,3.0e9,])
    Price_freq_range_norm = (Price_band/Price_central_freq)

    Dish_D = np.array([100,100,64])

    Price_BL_SEFD = np.array([calc_SEFD(calc_DishArea(Dish_D[0]), 20, eff=0.72),
        calc_SEFD(calc_DishArea(Dish_D[1]), 20, eff=0.72),
        calc_SEFD(calc_DishArea(Dish_D[2]), 35, eff=0.7),
        ])

    m=10.; nu =3.; t=300.
    Price_Sens = np.array([calc_Sensitivity(m,nu,t,SEFD=Price_BL_SEFD[i]) for i in range(len(Price_BL_SEFD))])

    dist_std = (50*3.26156 * u.lyr.to('m'))
    Price_BL_EIRP = np.array([calc_EIRP_min(dist_std,Sen) for Sen in Price_Sens])
    Price_BL_rarity = Price_BL_stars*Price_freq_range_norm

    Price_BL_stars_tot = Price_BL_stars[:2].sum()
    Price_rarity_tot = Price_BL_rarity[:2].sum()
    Price_EIRP_tot = Price_BL_EIRP[:2].max()

    iband = 900e6
    Price_BL_speed = Price_BL_SEFD.mean()**2*nu/iband

    Price_BL_sky = Price_BL_stars * np.array([calc_BeamSize(Dish_D[i],Price_central_freq[i]) for i in range(len(Dish_D))])
    Price_BL_DFM = Price_BL_sky * Price_band / Price_Sens**(3/2.)

    if verbose:
        print 'SEFD (Price_BL):', Price_BL_SEFD
        print 'Sens (Price_BL):', Price_Sens
        print 'EIRP (Price_BL):', Price_BL_EIRP
        print 'BeamSize (Price_BL):', np.array([calc_BeamSize(Dish_D[i],Price_central_freq[i]) for i in range(len(Dish_D))])
        print 'Sky Coverage (Price_BL):', Price_BL_sky.sum()
        print 'CWTFM (Price_BL):',  zeta_AO *(Price_BL_EIRP) / (Price_BL_rarity)
        print 'CWTFM (Price_BL_tot):',  zeta_AO *(Price_EIRP_tot) / (Price_rarity_tot)

        print 'DFM (Price_BL):', Price_BL_DFM
        print '~o~'

    #---------------------------------------------------------------------------------
    #EIRP values in watts.
    P = np.array([1e10,1e12,1e14,1e16,1e18,1e20,1e23])

    #---------------------------
    # Luminosity function limit on putative transmitters.
    plt.plot(np.log10(P),np.log10(calc_NP_law3(P)),lw=20,color='gray',alpha=0.3)#,label=r'$\alpha$: %s'%alpha)

    plt.plot([17,17],[-11,4],'--',lw=5,color='black',alpha=0.5)#,label='Kardashev Type I')
    plt.plot([13,13],[-11,4],lw=5,color='black',alpha=0.5)#,label='AO Planetary Radar')
    #---------------------------

    alpha = 0.7
    markersize = 20
    fontsize = 20
    ticksize = fontsize - 2
    dot_size = markersize - 12

    plt.plot([np.log10(Price_BL_EIRP[2])],[np.log10(1./Price_BL_rarity[2])],'h', color = '#690182',markeredgecolor='w',markersize = markersize, label='Price (2019 - Parkes)')
    plt.plot([np.log10(Price_EIRP_tot)],[np.log10(1./Price_rarity_tot)],'h', color = '#440154',markeredgecolor='w',markersize = markersize,label='Price (2019 - GBT)')

    #-------
    plt.plot([np.log10(BL_EIRP)],[np.log10(1./BL_rarity)],'h', color = 'orange',markeredgecolor='#1c608e',markersize = markersize, label='Enriquez (2017)')

    plt.plot(np.log10(GM_EIRP),np.log10(1./GM_rarity),'o',color ='#303b6b',markeredgecolor='w',markersize = markersize,label='Gray&Mooley (2017)')

    plt.plot(np.log10(ATA_EIRP[0:2]),np.log10(1./ATA_rarity[0:2]),'^',color ='#a1c625',markeredgecolor='w',markersize = markersize,label='Harp (2016) a,b')
    plt.plot(np.log10(ATA_EIRP[2]),np.log10(1./ATA_rarity[2]),'s',color ='#a1c625',markeredgecolor='w',markersize = markersize,label='Harp (2016) c')
    plt.plot(np.log10(ATA_EIRP[3]),np.log10(1./ATA_rarity[3]),'h',color ='#a1c625',markeredgecolor='w',markersize = markersize,label='Harp (2016) d')
    plt.plot(np.log10(ATA_EIRP[0]),np.log10(1./ATA_rarity_tot),'^w',markersize = markersize,markeredgecolor='#a1c625',alpha=0.7,label='Harp (2016) All*')
    plt.plot(np.log10(ATA_EIRP[0:1]),np.log10(1./ATA_rarity[0:1]),'ow',markersize = dot_size, markeredgecolor='#a1c625')

    plt.plot([np.log10(Siemion_EIRP)],[np.log10(1./Siemion_rarity)],'>',color ='#26828e',markeredgecolor='w',markersize = markersize,label='Siemion (2013)')
    plt.plot([np.log10(Siemion_EIRP)],[np.log10(1./Siemion_rarity)],'ow',markersize = dot_size,markeredgecolor='#440154')

    plt.plot(np.log10(Ph_EIRP),np.log10(1./Ph_rarity),'<b',color ='#31688e',markeredgecolor='w',markersize = markersize,label='Phoenix')
    plt.plot([np.log10(Ph_EIRP_tot)],[np.log10(1./Ph_rarity_tot)],'<w',markeredgecolor='#31688e',markersize = markersize,label='Phoenix All*')

    plt.plot(np.log10(Horowitz_EIRP),np.log10(1./Horowitz_rarity),'oc',color ='#add8e6',markeredgecolor='w',markersize = markersize,label='Horowitz&Sagan (1993)')
    plt.plot(np.log10(Valdes_EIRP),np.log10(1./Valdes_rarity),'sy',color ='pink',markeredgecolor='w',markersize = markersize,label='Valdes (1986)')

    plt.plot([np.log10(Tarter_EIRP)],[np.log10(1./Tarter_rarity)],'vc',color ='#1f9e89',markeredgecolor='w',markersize = markersize,label='Tarter (1980)')
    plt.plot(np.log10(Verschuur_EIRP),np.log10(1./Verschuur_rarity),'sm',color ='#efda21',markeredgecolor='w',markersize = markersize,label='Verschuur (1973)')

    plt.xlabel('EIRP [log(W)]',fontsize = fontsize)
    #plt.ylabel('Transmiter Galactic Rarity [log((Nstars*BW)^-1)]',fontsize=fontsize)
    #plt.ylabel('Transmitter Rate \n [log(1/(Nstars * rel_BW))]',fontsize=fontsize)
    plt.ylabel('Transmitter Rate ',fontsize=fontsize)

    plt.xticks(fontsize = ticksize)
    plt.yticks(fontsize = ticksize)
    plt.ylim(-10,4)
    plt.xlim(10,23)

def compare_SETI_limits(EIRP,rarity,shape='o',color='k',project='This Project'):
    ''' Compare SETI project with previus surveys.
    '''

    #---------------------------------------------------------------------------------
    # plotting setup
    plt.ion()
    plt.figure(figsize=(15, 10))
    markersize = 20

    plt.plot([np.log10(EIRP)],[np.log10(1./rarity)],shape, color = color,markersize = markersize, label=project)
    ET_power_law()
    plt.legend(numpoints=1,scatterpoints=1,fancybox=True, shadow=True)

    plt.savefig('SETI_limits_comparison.png', format='png',bbox_inches='tight')
#     plt.savefig('Transmitter_Rarity_FoM.pdf', format='pdf', dpi=300,bbox_inches='tight')


