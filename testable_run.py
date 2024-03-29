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
from village import Village

from calc_outstats import allelefreq_fx
from plotting import (plot_allele_frequency,
        plot_coordinates_host)

from IPython import embed
import matplotlib.pyplot as plt



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
    # villages = [Village(**i) for i in ]

    # simulation
    sh = "simulation"
    burn_in = config.getint(sh, "burn_in")

    # host_demography
    sh = 'host_demography'
    villages = config.getint(sh, 'villages')
    hostpopsize = list(map(int, config.get(sh, 'hostpopsize').split(",")))
    prevalence = list(map(float, config.get(sh, 'prevalence').split(",")))
    muWormBurden = list(map(int, config.get(sh, 'muWormBurden').split(",")))
    sizeWormBurden = list(
        map(int, config.get(sh, 'sizeWormBurden').split(",")))
    assert villages == len(hostpopsize)
    muTrans = config.getint(sh, 'muTrans')
    sizeTrans = config.getint(sh, 'sizeTrans')

    # host_demography - between village parameters
    hostmigrate = config.getfloat(sh, 'hostmigrate')
    initial_distance_m = list(
        map(int, config.get(sh, 'initial_distance_m').split(",")))
    assert len(initial_distance_m) == (villages * (villages - 1)) / 2

    # vector
    sh = 'vector'
    sigma = config.getfloat(sh, 'sigma')
    bitesPperson = list(map(int, config.get(sh, 'bitesPperson').split(",")))
    hours2bite = list(map(int, config.get(sh, 'hours2bite').split(",")))
    densitydep_uptake = config.getboolean(sh, "densitydep_uptake")

    # parasite
    sh = 'parasite'
    fecund = config.getint(sh, 'fecund')
    surv_Juv = config.getfloat(sh, 'surv_Juv')
    shapeMF = config.getfloat(sh, 'shapeMF')
    scaleMF = config.getfloat(sh, 'scaleMF')
    shapeAdult = config.getfloat(sh, 'shapeAdult')
    scaleAdult = config.getfloat(sh, 'scaleAdult')
    densitydep_surv = config.getboolean(sh, 'densitydep_surv')
    densitydep_fec = config.getboolean(sh, 'densitydep_fec')

    # genetic
    sh = 'genetic'
    locus = config.getint(sh, 'locus')
    initial_migration = config.getfloat(sh, 'initial_migration')
    theta = [list(map(int,i.split(","))) for i in config.get(sh, 'theta').split()]
    basepairs = list(map(int, config.get(sh, 'basepairs').split(",")))
    mutation_rate = list(map(float, config.get(sh, 'mutation_rate').split(",")))
    recombination_rate = list(map(float, config.get(sh, 'recombination_rate').split(",")))
    time2Ancestral = config.getint(sh, 'time2Ancestral')
    thetaRegional = config.getfloat(sh, 'thetaRegional')
    time_join = config.getint(sh, 'time_join')
    selection = config.getboolean(sh, 'selection')
    fitness = config.getint(sh, 'fitness')
    perc_locus = list(map(float, config.get(sh, 'perc_locus').split(",")))
    cds_length = config.getint(sh, 'cds_length')
    intgen_length = config.getint(sh, 'intgen_length')

    # treatment
    sh = 'treatment'
    bednets = config.getboolean(sh, 'bednets')
    bnstart = list(map(int, config.get(sh, 'bnstart').split(",")))
    bnstop = list(map(int, config.get(sh, 'bnstop').split(",")))
    bncoverage = list(map(float, config.get(sh, 'bncoverage').split(",")))
    mda = config.getboolean(sh, 'mda')
    mda_start = config.getint(sh, 'mda_start')
    mda_num = config.getint(sh, 'mda_num')
    mda_freq = config.getint(sh, 'mda_freq')
    mda_coverage = list(map(float, config.get(sh, 'mda_coverage').split(",")))
    mda_macro = config.getfloat(sh, 'mda_macro')
    mda_micro = config.getfloat(sh, 'mda_micro')
    mda_sterile = config.getfloat(sh, 'mda_sterile')
    mda_clear = config.getint(sh, 'mda_clear')

    # outfiles
    sh = 'outfiles'
    demofigs=config.getboolean(sh, 'demofigs')
    demotables=config.getboolean(sh, 'demotables')
    demoTime=config.getint(sh, 'demoTime')
    sample_size=config.getfloat(sh, 'sample_size')
    window_length=config.getint(sh, 'window_length')
    num_windows=config.getint(sh, 'num_windows')
    popgenfigs=config.getboolean(sh, 'popgenfigs')
    popgentables=config.getboolean(sh, 'popgentables')
    popgenTime=config.getint(sh, 'popgenTime')
    wb2vcf=config.getboolean(sh, 'wb2vcf')
    wbmfvcf=config.getboolean(sh, 'wbmfvcf')
    wbadultvcf=config.getboolean(sh, 'wbadultvcf')
    wbjuvvcf=config.getboolean(sh, 'wbjuvvcf')
    wbfracvcf=config.getfloat(sh, 'wbfracvcf')
    wb2scikit=config.getboolean(sh, 'wb2scikit')
    wbdebug=config.getboolean(sh, 'wbdebug')

    if selection:
        if mda:
            if fitness == 1:
                 from survival_mda import survivalmda_sel1_fx as survfx
            elif fitness == 2:
                 from survival_mda import survivalmda_sel2_fx as survfx
            else:
                 from survival_mda import survivalmda_fx as survfx
        else:
            from survival import survivalbase_fx as survfx
    else:
        if mda:
            from survival_mda import survivalmda_fx as survfx
        else:
            from survival import survivalbase_fx as survfx
    mdalist = [mda_start, mda_num, mda_freq, mda_coverage, mda_macro, mda_micro, mda_sterile, mda_clear]
    bnlist = [bednets, bnstart, bnstop, bncoverage]
    cdslist = [perc_locus, cds_length, intgen_length]
    # set counters
    sim_time = numberGens

    dfHost, dfAdult, dfJuv, dfMF, dfSel, cds_coordinates =\
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
                               cdslist)

    ## after intialize run main loop
    fig = plot_coordinates_host(dfHost)
    embed()
    for month in range(sim_time):
        print("month is %i\n\n" %month)
        dfHost, dfJuv, dfMF, L3trans = trans.transmission_fx(month,
                                                            villages,
                                                            hostpopsize,
                                                            sigma,
                                                            bitesPperson,
                                                            hours2bite,
                                                            densitydep_uptake,
                                                            bnlist,
                                                            dfHost,
                                                            dfJuv,
                                                            dfMF)
        dfHost, dfAdult, dfJuv, dfMF, dfSel = survfx(month,
                                                    surv_Juv,
                                                    shapeMF,
                                                    scaleMF,
                                                    shapeAdult,
                                                    scaleAdult,
                                                    fecund,
                                                    locus,
                                                    mutation_rate,
                                                    recombination_rate,
                                                    basepairs,
                                                    selection,
                                                    dfSel,
                                                    cds_coordinates,
                                                    hostmigrate,
                                                    mdalist,
                                                    densitydep_surv,
                                                    densitydep_fec,
                                                    dfHost,
                                                    dfAdult,
                                                    dfJuv,
                                                    dfMF)
        if month > burn_in:
                  allelefreq_fx(dfAdult, dfSel)
                  dfAdult.groupby("village").describe()
                  dfJuv.groupby("village").describe()
                  dfMF.groupby("village").describe()
                  dfHost.groupby("village").describe()

    return(dfHost,dfAdult,dfJuv,dfMF,dfSel)

if __name__ == '__main__':
     # this probably needs to be run for at least 240 - 360 months to get away from starting conditions
     dfHost, dfAdult, dfJuv, dfMF, dfSel = wb_sims(10, 'tests/wbsims.cfg')

