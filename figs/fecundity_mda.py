#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
from figs.recombination import recombination_fx
from figs.mutation import mutation_fx
from figs.selection import selection_fx

def fecunditymda_fx(villages,
                    fecund,
                    locus,
                    mutation_rate,
                    recombination_rate,
                    basepairs,
                    selection,
                    cds_coordinates,
                    densitydep_fec,
                    clear_count,
                    mda_sterile,
                    mda_clear,
                    dfHost,
                    dfAdult,
                    dfSel):

    '''function for reduced fecundity under mda
    conditions: mda=True, selection=False

    Parameters
    ----------
    fecund: int
         rate of the poisson distribution for offspring number
    sterile_p: float
         proportion of the adult parasite that are permanently sterilized from drugs
    clear_time: int
         number of months the drugs is effective or present in the host
    clear_count: int
         count of months since drug administered
    dfAdult: dataframe
         pandas dataframe containing adult parasites
    dfMF: dataframe
          pandas dataframe containing larval stage parasites

    Returns
    ------
    dfAdult
    dfMF
    '''
    #by host or by village?
    if clear_count == 1: #permanent sterility
         for index, row in dfHost[dfHost.MDA == 1].iterrows():
             try:
                 #randomly select dfAdults, change sex to "S" for sterile
                 sterile = dfAdult.loc[dfAdult.hostidx == row.hostidx].sample(frac = mda_sterile)["sex"].index
                 dfAdult.ix[sterile,"sex"] = "S"
             except ValueError:
                  pass
    if clear_count > 0 and clear_count <= mda_clear: #Drugs cause temporary sterility over clear_time
        #linear function defining fecundity during drug clearance
        mmda = float(fecund - 1) / (mda_clear - 1 )
        bmda = 1 - mmda * 1
        #new base fecundity under drugs
        sterile_t = (mmda * clear_count + bmda)
        #assign value to dfAdult.fec
        dfAdult.loc[dfAdult.age < 6, "fec"] = np.random.poisson(sterile_t,
                    len(dfAdult[dfAdult.age < 6]))

        #linear function defining decline in fecundity with age
        m = float(0 - sterile_t) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        dfAdult.loc[dfAdult.age >= 6, "fec"] = np.random.poisson(m
                  * dfAdult.loc[dfAdult.age >= 6,"age"] + b, len(dfAdult[dfAdult.age >= 6]))
    else: #base fecundity when no drugs
        #all locations where age is less than 6
        dfAdult.loc[dfAdult.age < 6, "fec"] = np.random.poisson(fecund,
                  len(dfAdult[dfAdult.age < 6]))
        #linear function defining decline in fecundity with age
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        dfAdult.loc[dfAdult.age >= 6, "fec"] = np.random.poisson(m
                  * dfAdult.loc[dfAdult.age >= 6,"age"] + b, len(dfAdult[dfAdult.age >= 6]))

    #sex, recombination, mutation
    dfAdult_mf = recombination_fx(locus, dfAdult, recombination_rate, basepairs)
    dfAdult_mf, positions = mutation_fx(locus, dfAdult_mf, mutation_rate, recombination_rate, basepairs)
    if selection:
        dfAdult_mf, dfSel = selection_fx(dfAdult_mf, positions, dfSel, locus, cds_coordinates)

    return(dfAdult_mf, dfSel)

def fecunditymda_sel1_fx(villages,
                    fecund,
                    locus,
                    mutation_rate,
                    recombination_rate,
                    basepairs,
                    selection,
                    dfSel,
                    cds_coordinates,
                    densitydep_fec,
                    clear_count,
                    mda_sterile,
                    mda_clear,
                    dfHost,
                    dfAdult):

    '''function for reduced fecundity under mda option 1
    option 1 simplifies that when no MDA or selective event all phenotypes
    are essetially wildtype, so fitness is not evaluated
    conditions: mda=True, selection=True, 1

    Parameters
    ----------
    fecund: int
         rate of the poisson distribution for offspring number
    sterile_p: float
         proportion of the adult parasite that are permanently sterilized from drugs
    clear_time: int
         number of months the drugs is effective or present in the host
    clear_count: int
         count of months since drug administered
    dfAdult: dataframe
         pandas dataframe containing adult parasites
    dfMF: dataframe
          pandas dataframe containing larval stage parasites

    Returns
    ------
    dfAdult
    dfMF
    '''
    if clear_count == 1: #permanent sterility
         for index, row in dfHost[dfHost.MDA == 1].iterrows():
             #randomly select dfAdults, change sex to "S" for sterile
             try:
                 mdarand = np.random.random(len(dfAdult.hostidx == row.hostidx))
                 mdasterile = dfHost.loc[np.where(mdarand < mda_sterile
                              ** dfAdult[dfAdult.hostidx == row.hostidx]["selF"])].index
                 dfAdult.ix[mdasterile, "sex"] = "S"
             except ValueError:
                 pass
    if clear_count > 0 and clear_count <= mda_clear: #Drugs cause temporary sterility over clear_time
         #linear function defining fecundity during drug clearance
         mmda = float(fecund - 1) / (mda_clear - 1 )
         bmda = 1 - mmda * 1
         #new base fecundity under drugs
         sterile_t = (mmda * clear_count + bmda)
         #assign value to dfAdult.fec
         dfAdult.loc[dfAdult.age < 6, "fec"] = np.random.poisson(sterile_t **
                    dfAdult.loc[dfAdult.age < 6, "selF"],
                               len(dfAdult[dfAdult.age < 6]))

        #linear function defining decline in fecundity with age
         m = float(0 - sterile_t) / (21 - 6)
         b = 0 - m * 21
         #assign fecundity value based on age function
         dfAdult.loc[dfAdult.age >= 6, "fec"] = np.random.poisson((m
              * dfAdult.loc[dfAdult.age >= 6,"age"] + b) ** dfAdult.loc[dfAdult.age >= 6,"selF"],
                    len(dfAdult[dfAdult.age >= 6]))

    else: #base fecundity when no drugs
        #all locations where age is less than 6
        dfAdult.loc[dfAdult.age < 6, "fec"] = np.random.poisson(fecund,
                  len(dfAdult[dfAdult.age < 6]))
        #linear function defining decline in fecundity with age
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        dfAdult.loc[dfAdult.age >= 6, "fec"] = np.random.poisson(m
                  * dfAdult.loc[dfAdult.age >= 6,"age"] + b, len(dfAdult[dfAdult.age >= 6]))

    #sex, recombination, mutation
    dfAdult_mf = recombination_fx(locus, dfAdult, recombination_rate, basepairs)
    dfAdult_mf, positions = mutation_fx(locus, dfAdult_mf, mutation_rate, recombination_rate, basepairs)
    if selection:
        dfAdult_mf, dfSel = selection_fx(dfAdult_mf, positions, dfSel, locus, cds_coordinates)

    return(dfAdult_mf, dfSel)

def fecunditymda_sel2_fx(villages,
                    fecund,
                    locus,
                    mutation_rate,
                    recombination_rate,
                    basepairs,
                    selection,
                    dfSel,
                    cds_coordinates,
                    densitydep_fec,
                    clear_count,
                    mda_sterile,
                    mda_clear,
                    dfHost,
                    dfAdult):
    '''function for reduced fecundity under mda option 2
    option 2 is when the mutant are less fit then the wildtype when no mda
    is being applied.
    conditions: mda=True, selection=True, 2

    Parameters
    ----------
    fecund: int
         rate of the poisson distribution for offspring number
    sterile_p: float
         proportion of the adult parasite that are permanently sterilized from drugs
    clear_time: int
         number of months the drugs is effective or present in the host
    clear_count: int
         count of months since drug administered
    dfAdult: dataframe
         pandas dataframe containing adult parasites
    dfMF: dataframe
          pandas dataframe containing larval stage parasites

    Returns
    ------
    dfAdult
    dfMF
    '''
    if clear_count == 1: #permanent sterility
         for index, row in dfHost[dfHost.MDA == 1].iterrows():
             #randomly select dfAdults, change sex to "S" for sterile
             try:
                 mdarand = np.random.random(len(dfAdult.hostidx == row.hostidx))
                 mdasterile = dfHost.loc[np.where(mdarand < mda_sterile
                              ** dfAdult[dfAdult.hostidx == row.hostidx]["selF"])].index
                 dfAdult.ix[mdasterile, "sex"] = "S"
             except ValueError:
                 pass

    if clear_count > 0 and clear_count <= mda_clear: #Drugs cause temporary sterility over clear_time
         #linear function defining fecundity during drug clearance
         mmda = float(fecund - 1) / (mda_clear - 1 )
         bmda = 1 - mmda * 1
         #new base fecundity under drugs
         sterile_t = (mmda * clear_count + bmda)
         #assign value to dfAdult.fec
         dfAdult.loc[dfAdult.age < 6, "fec"] = np.random.poisson(sterile_t **
                    dfAdult.loc[dfAdult.age < 6, "selF"],
                               len(dfAdult[dfAdult.age < 6]))

        #linear function defining decline in fecundity with age
         m = float(0 - sterile_t) / (21 - 6)
         b = 0 - m * 21
         #assign fecundity value based on age function
         dfAdult.loc[dfAdult.age >= 6, "fec"] = np.random.poisson((m
              * dfAdult.loc[dfAdult.age >= 6,"age"] + b) ** dfAdult.loc[dfAdult.age >= 6,"selF"],
                    len(dfAdult[dfAdult.age >= 6]))
    else: #base fecundity when no drugs
         #linear function defining fecundity during drug clearance
         mmda = float(fecund - 1) / (mda_clear - 1 )
         bmda = 1 - mmda * 1
         #new base fecundity under drugs
         sterile_t = (mmda * clear_count + bmda)
         #assign value to dfAdult.fec
         dfAdult.loc[dfAdult.age < 6, "fec"] = np.random.poisson(sterile_t **
                    (1 - abs(1 - dfAdult.loc[dfAdult.age < 6, "selF"])),len(dfAdult[dfAdult.age < 6]))

         #linear function defining decline in fecundity with age
         m = float(0 - sterile_t) / (21 - 6)
         b = 0 - m * 21
         #assign fecundity value based on age function
         dfAdult.loc[dfAdult.age >= 6, "fec"] = np.random.poisson((m
              * dfAdult.loc[dfAdult.age >= 6,"age"] + b) ** (1 - abs(1 - dfAdult.loc[dfAdult.age >= 6,"selF"])),
                 len(dfAdult[dfAdult.age >= 6]))

    #sex, recombination, mutation
    dfAdult_mf = recombination_fx(locus, dfAdult, recombination_rate, basepairs)
    dfAdult_mf, positions = mutation_fx(locus, dfAdult_mf, mutation_rate, recombination_rate, basepairs)
    if selection:
        dfAdult_mf, dfSel = selection_fx(dfAdult_mf, positions, dfSel, locus, cds_coordinates)

    return(dfAdult_mf, dfSel)
