# -*- coding: utf-8 -*-
"""
Spyder Editor

"""
import numpy as np
import pandas as pd
import math
from sklearn.metrics import pairwise_distances
import random
from scipy.stats import weibull_min
import subprocess

import wbsims_initialize as wbinit
import transmission as trans
from survival import survivalbase_fx
   
 
def wb_sims(numberGens):
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
    #simulation
    numberGens = 1000
    burn_in = 360
    #host_demography
    villages = 2
    hostpopsize = [100, 200]
    prevalence = [0.1, 0.3]
    hostmigrate = 0
    muTrans = 100
    sizeTrans = 1
    muWormBurden = [5, 5]
    sizeWormBurden = [50, 50]
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

    #set counters
    month = 0
    sim_time = numberGens
    
    if selection:
         dfAdult, dfHost, dfMF, dfJuv, dfSel = wbinit.wbsims_init(villages, hostpopsize, prevalence, muTrans, sizeTrans, 
                                  muWormBurden, sizeWormBurden, locus, initial_migration, 
                                  initial_distance_m, theta, basepairs, mutation_rate, 
                                  recombination_rate, time2Ancestral, thetaRegional,
                                  time_join, selection)
         while month <= sim_time:
              dfJuv, dfHost = trans.transmission_fx(villages, hostpopsize, sigma, bitesPperson, 
                                        hours2bite, densityDep, bednets, bnstart,
                                        bnstop, bncoverage, month, dfMF, dfJuv, dfHost)
              dfAdult, dfJuv, dfMF, dfHost, dfSel = survivalbase_fx(month, surv_Juv, shapeMF, scaleMF, shapeAdult,
                                                   scaleAdult, dfMF, dfAdult, dfJuv, dfHost,
                                                   fecund, locus, mutation_rate, recombination_rate, 
                                                   basepairs, selection, dfSel) 
              month += 1

         
    else:
         dfAdult, dfHost, dfMF, dfJuv = wbinit.wbsims_init(villages, hostpopsize, prevalence, muTrans, sizeTrans, 
                                  muWormBurden, sizeWormBurden, locus, initial_migration, 
                                  initial_distance_m, theta, basepairs, mutation_rate, 
                                  recombination_rate, time2Ancestral, thetaRegional,
                                  time_join, selection)        
         while month <= sim_time:
             dfJuv, dfHost = trans.transmission_fx(villages, hostpopsize, sigma, bitesPperson, 
                                             hours2bite, densityDep, bednets, bnstart,
                                             bnstop, bncoverage, month, dfMF, dfJuv, dfHost)
             dfAdult, dfJuv, dfMF, dfHost = survivalbase_fx(month, surv_Juv, shapeMF, scaleMF, shapeAdult,
                                                        scaleAdult, dfMF, dfAdult, dfJuv, dfHost,
                                                        fecund, locus, mutation_rate, recombination_rate, 
                                                        basepairs, selection, hostmigrate) 
             month += 1

if __name__ == '__main__':
     #this probably needs to be run for at least 240 - 360 months to get away from starting conditions
     wb_sims(24)
 
