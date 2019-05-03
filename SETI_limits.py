import numpy as np
from astropy import units as u
from scipy import stats
import pandas as pd

import pdb



class Telescope:

    def __init__(self,telescope_name='',diameter=None,Tsys=None,eff=None,station_area=None,N_stations=1):
        """ If array => diameter == baseline
        """

        self.telescope_name = telescope_name
        self.diameter = diameter
        self.Tsys = Tsys
        self.eff = eff
        self.station_area = station_area
        self.N_stations = N_stations

        if self.N_stations == 1 and self.diameter:
            self.area_eff = self.calc_DishArea()
        elif self.station_area and self.N_stations > 1:
            self.area_eff = self.station_area * self.N_stations
        else:
            print 'Need dish diameter or total area.'

    def calc_DishArea(self):
        """ Compute dish area
        d = dish diameter
        """

        return np.pi * (self.diameter/2.)**2

    def calc_gain(self,ll, d):
        """ Gain of a dish telescope
        ll = wavelength (lambda)
        d = dish diameter (m)
        """

        return (np.pi * d / ll)**2

    def calc_Tsky_low_freq(self,v):
        ''' Tsky calculation at low frequencies ().
        '''

        c = 2.998e8  #speed of light
        l = c / v
        T0 = 60   # in K
        Tsky = T0 * l**2.55

        return Tsky


class SETI_project:

    def __init__(self, telescope=None, N_stars=None, N_pointings=None, band = None, inst_band = None, central_freq = None, snr = None, delta_freq = None, obs_time = None, distance = None):

        self.telescope = telescope
        self.N_stars = N_stars
        if N_pointings:
            self.N_pointings = N_pointings
        else:
            self.N_pointings = N_stars

        self.band = band
        self.inst_band = inst_band
        self.central_freq = central_freq

        self.snr = snr
        self.delta_freq = delta_freq
        self.obs_time = obs_time

        self.distance = self.init_distance(distance)

        self.N_stars_perBeam = None


    def calc_standard_metric(self,):
        #---------------------------
        # Standardizing the sensitivity of telescopes by figuring out the max EIRP of a transmitter they could have found.
        # For this, I need the sensitivity of the observation, given SEFD and so on.

        #Standard distance (100 ly)
        dist_std = (100. * u.lyr.to('m'))
        #So decided not to use the distance above, it does makes sense if the original distance is shorter, but not if the distance is farther away (like in Siemion2013)


        #nomalization to L_AO = 2e13 W, fraq_freq = 1/2 and N_stars = 1k
        zeta_AO  =  1 #np.log10(1e3*.5)/ np.log10(2e13)
        zeta_AO  =  1e3*0.5/ 1e13


    def init_distance(self,distance,pc=True):
        '''Init distance from light years or pc
        '''

        if not distance:
            dist_std = None
        else:
            #Set from either pc or ly.
            if pc:
                dist_std = (distance*3.26156 * u.lyr.to('m'))
            else:
                dist_std = (distance * u.lyr.to('m'))

        return dist_std

    def calc_SEFD(self,A, Tsys, eff=1.0):
        """ Calculate SEFD
        Tsys = system temperature
        A = collecting area
        Ae = effective collecting area
        eff = aperture efficency (0.0 to 1.0)
        """

        kb = 1.3806488e3  # 1.38064852e-23    Boltzmann constant
        Ae = A*eff

        sefd = 2 * Tsys * kb / Ae

        #Calculating factor for total SEFD of array, from SEFD of single station.
        if self.telescope.N_stations > 2:
            array_sefd_factor = self.telescope.N_stations * (self.telescope.N_stations-1)
        else:
            array_sefd_factor = 1.

        sefd =  sefd / np.sqrt(array_sefd_factor)

        return sefd

    def calc_Sensitivity(self,m, nu, t, sefd=None, Tsys=10,eff=1.0,A=1,npol=2.,nu_t=1.,narrow=True):
        """ Minimum detectable luminosity for narrowband emission
        Tsys = system temperature
        A = collecting area
        m = threshold, (e.g. 10)
        nu = channel bandwidth [Hz]
        t = observing time [sec]
        narrow = True if signal is narrower than spectral resolution.
        """

        if not sefd:
            SEFD = self.calc_SEFD(self.telescope.station_area,self.telescope.Tsys,eff=self.telescope.eff)
        else:
            SEFD = sefd

        if narrow:
            sens = m * SEFD * np.sqrt(nu/(npol*t)) / nu_t
        else:
            sens = m * SEFD / np.sqrt(npol*nu*t)

        return sens

    def calc_BeamSize(self,d,v,verbose=False):
        """ Compute BeamSize
        d = dish diameter
        v = frequency
        """

        c = 2.998e8  #speed of light [m/s]
        l = c/v

        #This factor is used for LOFAR HBA.
        #For the case of GBT, this is the first order, there is a second term which is frequency dependent, we are ignoring that here.
        factor = 1.02

        if 'LOFAR' in self.telescope.telescope_name and self.central_freq < 1e8:
            factor = 1.1    #Different factor used for the  LBAs

       #This is the old equation I was using.
       #The 1.22 came from wikipidea... which is more use for optical telescopes.
       #BeamSize = (1.22* (c/(d*v)) *57.2958)

        BeamSize = factor*l/d*57.2958

        if verbose:
            print '\nBeam is: %f \n'%(BeamSize)

        return BeamSize

    def calc_min_EIRP(self,d,Sens=None,verbose=False):
        """ Minimum detectable luminosity (EIRP) for narrowband emission.
        d = distance to target star []
        Sens = sensitivity of the obs (Jy)
        """

        unit_Jy2W = 1e-26 #1 Jy = 1e-26 W/m2/Hz)

        if not Sens:
            SEFD = self.calc_SEFD(self.telescope.station_area,self.telescope.Tsys,eff=self.telescope.eff)
            Sens = self.calc_Sensitivity(self.snr, self.delta_freq,self.obs_time,sefd=SEFD)

        EIRP  =  4 * np.pi * d**2 * Sens * unit_Jy2W

        if verbose:
            print 'Sens = ', Sens
            print 'SEFD = ', SEFD
            print 'EIRP = ', EIRP

        return EIRP

    def calc_sky_coverage(self):

        fwhm = self.calc_BeamSize(self.telescope.diameter,self.central_freq)
        sky = self.N_pointings * (fwhm/2.)**2*np.pi

        return sky

    def calc_DFM(self,Sky=None):
        ''' Calculate the Drake Figure of Merit
        '''

        SEFD = self.calc_SEFD(self.telescope.station_area,self.telescope.Tsys,eff=self.telescope.eff)
        Sens = self.calc_Sensitivity(self.snr, self.delta_freq,self.obs_time,sefd=SEFD)

        if not Sky:
            Sky = self.calc_sky_coverage()

        unit_Jy2W = 1e-26 #1 Jy = 1e-26 W/m2/Hz)

        #Getting units of (Gz m**3 / W**3/2)
        Sky /= 41253.   #Over full sky in deg2
        band = self.band / 1e9
        Sens = Sens * unit_Jy2W

        DFM = Sky * band / Sens**(3/2.)

        return DFM

    def stars_in_beam(self,distance,verbose=False,zeta_AO=False):
        ''' Estimates the number of stars in a beam to a given distance.
        '''

        #Standard distance (10 pc)
        #dist_std = (10 *3.26156 * u.lyr.to('m'))   # this will replace the 10 below. once the distance given is taken to the correct units.

        if zeta_AO:
            fwhm = 3/60.  #Based on Arecibo at central frequency of Lband.
        else:
            fwhm = self.calc_BeamSize(self.telescope.diameter,self.central_freq)
        N_stars = 317 * (distance / 10. )**3 * ((fwhm/2.)**2*np.pi / 41253.)

        if verbose:
            print "The number of starts to a distance of %.2f pc, is a %.2f ' beam is %.2f ."%(distance,fwhm*60.,N_stars)

        self.N_stars_perBeam = N_stars

        return N_stars

    def calc_CWTFM(self,d_beam = None, verbose=True):
        ''' Calculate the Continuous Wave Transmiter Figure of Merit (CWTFM)
        '''

        #---------------------------
        #nomalization to L_AO = 2e13 W, fraq_freq = 1/2 and N_stars = 1k
        #zeta_AO  =  np.log10(1e3*.5)/ np.log10(2e13)
        #Also tried a few other combination of log and ln
        #Interestingly sqrt(EIRP) is distance dependent.
        if d_beam:
#            zeta_AO  =  1e3*self.stars_in_beam(d_beam,zeta_AO=True)*0.5/ np.log10(2e13)

            zeta_AO  =  1e3*self.stars_in_beam(d_beam,zeta_AO=True)*0.5/ np.sqrt(2e13)
        else:
#            zeta_AO  =  1e3*0.5/ np.log10(2e13)
            zeta_AO  =  1e3*0.5/ np.sqrt(2e13)

        #---------------------------

        SEFD = self.calc_SEFD(self.telescope.station_area,self.telescope.Tsys,eff=self.telescope.eff)
        Sens = self.calc_Sensitivity(self.snr, self.delta_freq,self.obs_time,sefd=SEFD)

        freq_range_norm = self.band/self.central_freq

        if d_beam:
            self.distance = self.init_distance(d_beam)
            EIRP = self.calc_min_EIRP(self.distance,Sens)
            rarity = self.stars_in_beam(d_beam)* self.N_pointings * freq_range_norm

        else:
            EIRP = self.calc_min_EIRP(self.distance,Sens)
            rarity = self.N_stars * freq_range_norm

#        CWTFM  = zeta_AO * np.log10(EIRP) / rarity
        CWTFM  = zeta_AO * np.sqrt(EIRP) / rarity

        return CWTFM

    def calc_survey_speed(self):
        ''' Calculate the survey speed
        '''

        SEFD = self.calc_SEFD(self.telescope.station_area,self.telescope.Tsys,eff=self.telescope.eff)
        survey_speed = SEFD**2 * self.delta_freq /self.inst_band

        return survey_speed


def freq_game(f1,f2):
    ''' Playing around with what parameter of frequencies makes more sense.
    '''

    central_f = (f2+f1)/2.

    print '(f2-f1)/ central_f',(f2-f1)/ central_f
    print '(f2-f1)/log(central_f)',(f2-f1)/ np.log10(central_f)
    print '(f2-f1)/f1',(f2-f1)/ f1
    print 'f2/f1',f2/f1
    print 'log(f2/f1)',np.log10(f2/f1)


    # didn't like f2/f1 for the case comparing 1-10GHz and 2-11GHz ...since should be more similar than they turn out to be.



