#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import pandas as pd
from scipy.stats import weibull_min
import random

from figs.fecundity import fecunditybase_fx
from figs.host_migration import hostmigration_fx

def survivalbase_fx(month,
                    village,
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
                    hostmigrate,
                    mdalist,
                    densitydep_surv,
                    densitydep_fec,
                    dfHost,
                    dfAdult,
                    dfJuv,
                    dfMF):


    '''Base survival function
    Parameters
    ---------
    month: int
         time in months of the simulation
    surv_Juv: float
         survival prob of juvenille life stage
    shapeMF: float
         shape parameter for weibull distribution of MF
    scaleMF: int
         scale parameter for weibull distribution of MF
    shapeAdult: float
         shape parameter for weibull distribution of Adults
    scaleAdult: int
         scale parameter for weibull distribution of MF
    dfMF : df
        dataframe of MF
    dfAdult : df
        dataframe of adult worms
    dfJuv : df
        dataframe of juvenille worms
    dfHost : df
        dataframe of hosts
    basepairs : int, list
        length of loci
    hostmigrate : float
        rate of migration per year between villages
    selection : boolean
        T/F for selection
    dfSel : df
        dataframe of cds positions and fitness
    cds_coordinates : list
        list of coding seq coordinates
    mdalist : list
        list of mda parameters
    densitydep_surv : boolean
    densitydep_fec : boolean
    Returns
    -------
    dfMF
    dfAdult
    dfJuv
    dfHost
    dfSel

    '''
    #adult worms and hosts are only evaluated per year
    if month%12 == 0:
        #Adult survival is based on weibull cdf
        surv_adultrand = np.random.random(len(dfAdult.meta))
        try:
            surv_adultfxage = weibull_min.cdf(dfAdult.meta.age, shapeAdult,loc=0,scale=scaleAdult)
        except TypeError:
            surv_adultfxage = weibull_min.cdf(0, shapeAdult,loc=0,scale=scaleAdult)
        surviveAdult = np.where(surv_adultrand <= (1 - surv_adultfxage))
        dfAdult.meta = dfAdult.meta.iloc[surviveAdult]
        dfAdult.meta.age = dfAdult.meta.age + 1 #2 - 21

        ##host survival is from act table
        dfHost = dfHost[dfHost.age < dfHost.agedeath]
        #remove all worms with dead host.hostidx from all dataframes
        dfAdult.meta = dfAdult.meta.loc[dfAdult.meta["hostidx"].isin(dfHost.hostidx)]
        dfJuv.meta = dfJuv.meta.loc[dfJuv.meta["hostidx"].isin(dfHost.hostidx)]
        dfMF.meta = dfMF.meta.loc[dfMF.meta["hostidx"].isin(dfHost.hostidx)]
        #add 1 year to all ages of hosts
        dfHost.age = dfHost.age + 1
        if hostmigrate != 0:
            dfHost = hostmigration_fx(dfHost, hostmigrate)

    ##Juv is exponential 0.866; surv_Juv
    #dont include age 0 which just moved from transmission fx
    dfJuv.meta.age += 1
    surv_juvrand = np.random.random(len(dfJuv.meta))
    surviveJuv = np.where(surv_juvrand <= surv_Juv)
    dfJuv.meta = dfJuv.meta.iloc[surviveJuv]

    ##MF is weibull cdf
    surv_mfrand = np.random.random(len(dfMF.meta))
    try:
        surv_mffxage = weibull_min.cdf(dfMF.meta.age,shapeMF,loc=0,scale=scaleMF)
    except TypeError:
        surv_mffxage = weibull_min.cdf(0,shapeMF,loc=0,scale=scaleMF)
    surviveMF = np.where(surv_mfrand <= (1 - surv_mffxage))
    dfMF.meta = dfMF.meta.loc[surviveMF]
    dfMF.meta.age = dfMF.meta.age + 1 #2 - 12
    dfMF.meta = dfMF.meta[dfMF.meta.age < 13] #hard cutoff at 12 months

    ##move Juv age 13 to adult age 1
    #dfJuv_new = pd.DataFrame({})
    dfJuv_new = dfJuv.meta[dfJuv.meta.age > 12].copy()
    #reset age to adult
    dfJuv_new.age = 1
    #increase R0net for next gen
    dfJuv_new.R0net = dfJuv_new.R0net + 1
    #append to adults
    dfAdult.meta = pd.concat([dfAdult.meta, dfJuv_new], ignore_index=True)
    #remove Juv age 13 from dfJuv
    dfJuv.meta = dfJuv.meta[dfJuv.meta.age <= 12]

    ##call to fecundity fx to deepcopy adult to dfMF age 1
    #fecundity calls mutation/recombination
    dfAdult_mf = fecunditybase_fx(fecund, dfAdult, locus, mutation_rate,
                                         recombination_rate, basepairs, selection,
                                         densitydep_fec)
    dfAdult_mf.meta.age = 1
    dfAdult_mf.meta.fec = 0
    dfAdult_mf.meta.sex = [random.choice("MF") for i in range(len(dfAdult_mf.meta))]
    dfMF.meta = pd.concat([dfMF.meta, dfAdult_mf.meta], ignore_index=True)

    return(dfHost, dfAdult, dfJuv, dfMF)
