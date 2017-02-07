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
from village import Village


from calc_outstats import allelefreq_fx
from plotting import plot_allele_frequency



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
    # simulation
    burn_in = 360

    # host_demography
    sh = 'host_demography'
    villages = config.getint(sh, 'villages')
    hostpopsize = list(map(int, config.get(sh, 'hostpopsize').split(",")))
    prevalence = list(map(float, config.get(sh, 'prevalence').split(",")))
    muWormBurden = list(map(int, config.get(sh, 'muWormBurden').split(",")))
    sizeWormBurden = list(
        map(int, config.get(sh, 'sizeWormBurden').split(",")))
    assert villages == len(hostpopsize)

    # Between village parameters
    hostmigrate = config.getint(sh, 'hostmigrate')
    muTrans = config.getint(sh, 'muTrans')
    sizeTrans = config.getint(sh, 'sizeTrans')
    initial_distance_m = list(
        map(int, config.get(sh, 'initial_distance_m').split(",")))
    print(initial_distance_m)
    assert len(initial_distance_m) + 1 == villages

    # vector
    sh = 'vector'
    sigma = config.getint(sh, 'sigma')
    bitesPperson = list(map(int, config.get(sh, 'bitesPperson').split(",")))
    hours2bite = list(map(int, config.get(sh, 'hours2bite').split(",")))
    densityDep = list(map(bool, config.get(sh, 'densitydep_uptake').split(",")))

    # parasite
    sh = 'parasite'
    fecund = config.getint(sh, 'fecund')
    surv_Juv = config.getfloat(sh, 'surv_Juv')
    shapeMF = config.getfloat(sh, 'shapeMF')
    scaleMF = config.getint(sh, 'scaleMF')
    shapeAdult = config.getfloat(sh, 'shapeAdult')
    scaleAdult = config.getint(sh, 'scaleAdult')

    # genetic
    sh = 'genetic'
    locus = config.getint(sh, 'locus')
    initial_migration = config.getfloat(sh, 'initial_migration')
    theta = [[5, 5], [10, 10]]
    basepairs = list(map(int, config.get(sh, 'basepairs').split(",")))
    mutation_rate = list(map(float, config.get(sh, 'mutation_rate').split(",")))
    recombination_rate = list(map(float, config.get(sh, 'recombination_rate').split(",")))
    time2Ancestral = config.getint(sh, 'time2Ancestral')
    thetaRegional = config.getint(sh, 'thetaRegional')
    time_join = config.getint(sh, 'time_join')
    selection = True 
    perc_locus = [0, 0.18]
    #cds_length = config.getint(sh, 'cds_length')
    intgen_length = 2500 # treatment
    bnstart = [0, 0]
    bnstop = [0, 0]
    bncoverage = [0, 0]
    mda = [False, False]
    mda_start = [12, 12]
    mda_num = [6, 6]  # how many mdas
    mda_freq = 12  # every 12 months
    mda_coverage = [0.8, 0.7]
    mda_macro = 0.05
    mda_micro = 0.95
    mda_sterile = 0.35
    mda_clear = 6
    # output
    perc_locus = [0.2, 0.5]
    cds_length = [1000, 2000]
    intgen_length = 600

    # set counters
    month = 0
    sim_time = numberGens

    # SELECTION
    if selection:
        dfAdult, dfHost, dfSel, dfJuv, dfMF, cds_coordinates =\
                 wbinit.wbsims_init(villages,
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
        print("******INITIALIZED********")
        from IPython import embed
        test = allelefreq_fx(dfAdult, dfSel)
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

    # NOT SELECTION
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
     # this probably needs to be run for at least 240 - 360 months to get away from starting conditions
     wb_sims(10, 'tests/wbsims.cfg')
 
