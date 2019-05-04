from ET_power_law import *

import matplotlib.pylab as plt

#---------------------------
# Your new values ( e.g. Enriquez 2017)

project = 'This project'
telescope = 'GBT'
N_stars = 692 *100
band = 660e6 #Hz
central_freq = 1.5e9 #Hz

dish_diam = 100  #meters
dish_Tsys = 20  #K
dish_app_eff = 0.72  #Apperture Efficiency

SNR_threshold =  25 # sigma above the mean
spectral_resolution = 3. # Hz
scan_obs_time = 300 # sec
max_distance = 50  # ly

iband = 800e6  # Hz  Instantaneous Bandwidth

shape = '*'
color = 'red'

#---------------------------
#Calculating limits

zeta_AO  =  1e3*0.5/ 1e13

freq_range_norm = (band/central_freq)
SEFD = calc_SEFD(calc_DishArea(dish_diam), dish_Tsys, eff=dish_app_eff) # 10 Jy (GBT)
Sens = calc_Sensitivity(SNR_threshold, spectral_resolution,scan_obs_time,SEFD=SEFD)

dist_m = (max_distance*3.26156 * u.lyr.to('m'))

EIRP = calc_EIRP_min(dist_m,Sens)

survey_rarity = N_stars*freq_range_norm
survey_speed = SEFD**2*spectral_resolution/iband
survey_sky = N_stars * calc_BeamSize(dish_diam,central_freq)
survey_DFM = survey_sky * band / Sens**(3/2.)

def print_project():
    print '~o~', project ,' (', telescope,') ', '~o~'
    print 'SEFD :', SEFD
    print 'Sens :', Sens
    print 'EIRP :', EIRP
    print 'BeamSize :', calc_BeamSize(dish_diam,central_freq)
    print 'Sky Coverage :', survey_sky
    print 'CWTFM :',  zeta_AO *(EIRP) / (survey_rarity)
    print 'DFM :', survey_DFM


print_project()

#---------------------------
#Comparing SETI limits

compare_SETI_limits(EIRP,survey_rarity,shape=shape,color=color,project=project)

