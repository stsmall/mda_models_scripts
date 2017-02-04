# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
import subprocess
import math
import random
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from scipy.stats import weibull_min

import wbsims_initialize as wbinit
import transmission as trans
from survival import survivalbase_fx
<<<<<<< HEAD
from village import Village


from IPython import embed


=======
from calc_outstats import allelefreq_fx
>>>>>>> 57571acd4194c8604703d0635aad0608c65d176f
   
 
def wb_sims(numberGens, config_file):
    '''main function for simulations
    Parameters
    ---------
    numberGens : int
         how many months to run the simulation
    burn-in : int
         burn-in before recording data
    config_file : file
         file with options to be read by configparser()
    Returns
    -------
    '''
    config = configparser.ConfigParser()
    config.read(config_file)

    # villages = [Village(hostpopsize = 100, prevalence = 0.1)] 
    #simulation
#    numberGens = 1000
    burn_in = 360

    #host_demography
    sh = 'host_demography'
    villages = config.getint(sh, 'villages')
    hostpopsize = list(map(int, config.get(sh, 'hostpopsize').split(",")))
    prevalence = list(map(float, config.get(sh, 'prevalence').split(",")))
    muWormBurden = list(map(int, config.get(sh, 'muWormBurden').split(",")))
    sizeWormBurden = list(map(int, config.get(sh, 'sizeWormBurden').split(",")))
    assert villages == len(hostpopsize)

    # Between village parameters
    hostmigrate = config.getint(sh, 'hostmigrate')
    muTrans = config.getint(sh, 'muTrans')
    muTrans = 100
    sizeTrans = 1
    initial_distance_m = [1000]

    #vector
    sigma = 100
    bitesPperson = [10, 10]
    hours2bite = [8, 8]
    densityDep = [True, True]  
    #parasite
    fecund = 20
    surv_Juv = 0.866
    shapeMF = 3.3
    scaleMF = 10
    shapeAdult = 3.8
    scaleAdult = 8
    #genetic
    locus = 2
    initial_migration = 0.0001
    theta = [[5, 5], [10, 10]]
    basepairs = [13000, 200000]
    mutation_rate = [7.6E-8, 2.9E-9]
    recombination_rate = [0, 2.9E-9]
    time2Ancestral = 1800
    thetaRegional = 23
    time_join = 240
    selection = False
    perc_locus = [0, 0.18]
    cds_length = 1100
    intgen_length = 2500
    #treatment
    bednets = [False, False]
    bnstart = [0, 0]
    bnstop = [0, 0]
    bncoverage = [0, 0]
    mda = [False, False]
    mda_start = [12, 12]
    mda_num = [6, 6] #how many mdas
    mda_freq = 12 #every 12 months
    mda_coverage = [0.8, 0.7] 
    mda_macro = 0.05
    mda_micro = 0.95
    mda_sterile = 0.35
    mda_clear = 6
    #output
    perc_locus = [0.2, 0.5]
    cds_length = [1000, 2000]
    intgen_length = 600


    #set counters
    month = 0
    sim_time = numberGens
    
    ##SELECTION
    if selection:
         dfAdult, dfHost, dfSel, dfJuv, dfMF, cds_coordinates = wbinit.wbsims_init(villages, 
                                                                                   hostpopsize, 
                                                                                   prevalence, 
                                                                                   muTrans, 
                                                                                   sizeTrans, 
                                                                                   muWormBurden,
                                                                                   sizeWormBurden, 
                                                                                   locus, 
                                                                                   initial_migration, 
                                                                                   initial_distance_m, 
                                                                                   theta,
                                                                                   basepairs, 
                                                                                   mutation_rate, 
                                                                                   recombination_rate, 
                                                                                   time2Ancestral, 
                                                                                   thetaRegional,
                                                                                   time_join, 
                                                                                   selection, 
                                                                                   perc_locus, 
                                                                                   cds_length, 
                                                                                   intgen_length)
         for month in range(sim_time):
              dfJuv, dfHost = trans.transmission_fx(villages, 
                                                    hostpopsize, 
                                                    sigma, 
                                                    bitesPperson, 
                                                    hours2bite, 
                                                    densityDep, 
                                                    bednets, 
                                                    bnstart,
                                                    bnstop, 
                                                    bncoverage, 
                                                    month, 
                                                    dfMF, 
                                                    dfJuv, 
                                                    dfHost)
              dfAdult, dfJuv, dfMF, dfHost, dfSel = survivalbase_fx(month, 
                                                                    surv_Juv, 
                                                                    shapeMF, 
                                                                    scaleMF, 
                                                                    shapeAdult,
                                                                    scaleAdult, 
                                                                    dfMF, 
                                                                    dfAdult, 
                                                                    dfJuv, 
                                                                    dfHost,
                                                                    fecund, 
                                                                    locus, 
                                                                    mutation_rate, 
                                                                    recombination_rate, 
                                                                    basepairs, 
                                                                    selection, 
                                                                    dfSel) 
              if month > burn_in:
                  allelefreq_fx(dfAdult, dfSel)
                  dfAdult.groupby("village").describe()
                  dfJuv.groupby("village").describe()
                  dfMF.groupby("village").describe()
                  dfHost.groupby("village").describe()

    ##NOT SELECTION     
    else:
         dfAdult, dfHost, dfMF, dfJuv = wbinit.wbsims_init(villages, 
                                                           hostpopsize, 
                                                           prevalence, 
                                                           muTrans, 
                                                           sizeTrans, 
                                                           muWormBurden, 
                                                           sizeWormBurden, 
                                                           locus, 
                                                           initial_migration, 
                                                           initial_distance_m, 
                                                           theta, 
                                                           basepairs, 
                                                           mutation_rate, 
                                                           recombination_rate, 
                                                           time2Ancestral, 
                                                           thetaRegional,
                                                           time_join, 
                                                           selection, 
                                                           perc_locus, 
                                                           cds_length, 
                                                           intgen_length)        
         for month in range(sim_time):
             dfJuv, dfHost = trans.transmission_fx(villages, 
                                                   hostpopsize, 
                                                   sigma, 
                                                   bitesPperson, 
                                                   hours2bite, 
                                                   densityDep, 
                                                   bednets, 
                                                   bnstart,
                                                   bnstop, 
                                                   bncoverage, 
                                                   month, 
                                                   dfMF, 
                                                   dfJuv, 
                                                   dfHost)
             dfAdult, dfJuv, dfMF, dfHost = survivalbase_fx(month, 
                                                            surv_Juv, 
                                                            shapeMF, 
                                                            scaleMF, 
                                                            shapeAdult,
                                                            scaleAdult, 
                                                            dfMF, 
                                                            dfAdult, 
                                                            dfJuv, 
                                                            dfHost,
                                                            fecund, 
                                                            locus, 
                                                            mutation_rate, 
                                                            recombination_rate, 
                                                            basepairs, 
                                                            selection, 
                                                            hostmigrate) 
             if month > burn_in:
                 dfAdult.groupby("village").describe()
                 dfJuv.groupby("village").describe()
                 dfMF.groupby("village").describe()
                 dfHost.groupby("village").describe()
                 

if __name__ == '__main__':
     #this probably needs to be run for at least 240 - 360 months to get away from starting conditions
     wb_sims(1000, 'tests/wbsims.cfg')
 
