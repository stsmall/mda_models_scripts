#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import ipdb
import numpy as np
from .recombination import recombination_fx
from .mutation import mutation_fx
from .selection import selection_fx

def fecunditymda_fx(fecund,
                    dfworm,
                    locus,
                    mutation_rate,
                    recombination_rate,
                    basepairs,
                    selection,
                    densitydep_fec,
                    mda_sterile,
                    clear_count,
                    mda_clear,
                    hostmda):

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
####################################################
    #by host or by village?
    if clear_count == 1: #permanent sterility
        Adultmda = dfworm.meta[(dfworm.meta["stage"] == "A") &
                               (dfworm.meta.hostidx.isin(hostmda))].index.values
        #wormmda = dfworm.meta[dfworm.hostidx.isin(hostmda)].index.values
        #Adultmda = dfworm.meta.ix[wormmda, "stage" == "A"].index.values
        adultsterile = dfworm.meta.ix[Adultmda].groupby("hostidx").apply(lambda y: y.sample(frac = mda_sterile).index.values)
        steriix = np.concatenate([i for i in adultsterile.values]).ravel()
        dfworm.meta.ix[steriix,"sex"] = "S"

    else: pass

    if clear_count > 0 and clear_count <= mda_clear: #Drugs cause temporary sterility over clear_time
        adiix = dfworm.meta[dfworm.meta.stage == "A"].index.values
        young = adiix[np.where(dfworm.meta.ix[adiix].age < 6)]
        old = adiix[np.where(dfworm.meta.ix[adiix].age >= 6)]
        #linear function defining fecundity during drug clearance
        mmda = float(fecund - 1) / (mda_clear - 1 )
        bmda = 1 - mmda * 1
        #new base fecundity under drugs
        sterile_t = (mmda * clear_count + bmda)
        dfworm.meta.ix[young, 'fec'] = np.random.poisson(sterile_t, young.shape[0])

       #linear function defining decline in fecundity with age
        m = float(0 - sterile_t) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        positive_lambda = (dfworm.meta.ix[old].age.values * m) + b
        positive_lambda[positive_lambda < 0] = 0
        dfworm.meta.ix[old, 'fec'] = np.random.poisson(positive_lambda).astype(np.int64)
########################################################
    else: #base fecundity when no drugs
        adiix = dfworm.meta[dfworm.meta.stage == "A"].index.values
        young = adiix[np.where(dfworm.meta.ix[adiix].age < 6)]
        old = adiix[np.where(dfworm.meta.ix[adiix].age >= 6)]
        dfworm.meta.ix[young, 'fec'] = np.random.poisson(fecund, young.shape[0])
        #linear function defining decline in fecundity with age
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign fecundity value based on age function
        positive_lambda = (dfworm.meta.ix[old].age.values * m) + b
        positive_lambda[positive_lambda < 0] = 0
        dfworm.meta.ix[old, 'fec'] = np.random.poisson(positive_lambda).astype(np.int64)
    #sex, recombination, mutation
    dfAdult_mf = recombination_fx(locus, dfworm, adiix, recombination_rate, basepairs)
    # Positions is just the new positions
    dfAdult_mf, new_positions = mutation_fx(locus, dfAdult_mf,
         mutation_rate, recombination_rate, basepairs)
    if selection: #dfAdult.sel will be updated here to same length as dfAdult_mf.pos
        dfAdult_mf, dfworm = selection_fx(dfworm, dfAdult_mf, new_positions)
    return(dfAdult_mf, dfworm)

#def fecunditymda_sel1_fx(fecund,
#                    dfAdult,
#                    locus,
#                    mutation_rate,
#                    recombination_rate,
#                    basepairs,
#                    selection,
#                    densitydep_fec,
#                    mda_sterile,
#                    clear_count,
#                    mda_clear,
#                    dfHost):
#    '''function for reduced fecundity under mda option 1
#    option 1 simplifies that when no MDA or selective event all phenotypes
#    are essetially wildtype, so fitness is not evaluated
#    conditions: mda=True, selection=True, 1
#
#    Parameters
#    ----------
#    fecund: int
#         rate of the poisson distribution for offspring number
#    sterile_p: float
#         proportion of the adult parasite that are permanently sterilized from drugs
#    clear_time: int
#         number of months the drugs is effective or present in the host
#    clear_count: int
#         count of months since drug administered
#    dfAdult: dataframe
#         pandas dataframe containing adult parasites
#    dfMF: dataframe
#          pandas dataframe containing larval stage parasites
#
#    Returns
#    ------
#    dfAdult
#    dfMF
#    '''
#############################################################
#    if clear_count == 1: #permanent sterility
#         for index, row in dfHost[dfHost.MDA == 1].iterrows():
#             #randomly select dfAdults, change sex to "S" for sterile
#             try:
#                 mdarand = np.random.random(len(dfAdult.meta.hostidx == row.hostidx))
#                 mdasterile = dfHost.meta.ix[np.where(mdarand < mda_sterile
#                              ** dfAdult.meta[dfAdult.meta.hostidx == row.hostidx]["selF"])].index.values
#                 dfAdult.meta.ix[mdasterile, "sex"] = "S"
#             except ValueError:
#                 pass
#
#    if clear_count > 0 and clear_count <= mda_clear: #Drugs cause temporary sterility over clear_time
#         #linear function defining fecundity during drug clearance
#         mmda = float(fecund - 1) / (mda_clear - 1 )
#         bmda = 1 - mmda * 1
#         #new base fecundity under drugs
#         sterile_t = (mmda * clear_count + bmda)
#         assert sterile_t > 0
#         #assign value to dfAdult.fec
#         dfAdult.meta.loc[dfAdult.meta.age < 6, "fec"] = np.random.poisson(sterile_t **
#                    dfAdult.meta.loc[dfAdult.meta.age < 6, "selF"],
#                               len(dfAdult.meta[dfAdult.meta.age < 6]))
#
#        #linear function defining decline in fecundity with age
#         m = float(0 - sterile_t) / (21 - 6)
#         b = 0 - m * 21
#         #assign fecundity value based on age function
#         positive_lambda = (dfAdult.meta.loc[dfAdult.meta.age >= 6, "age"].values * m) + b
#         positive_lambda[positive_lambda < 0] = 0
#         dfAdult.meta.loc[dfAdult.meta.age >= 6, "fec"] = np.random.poisson((positive_lambda)
#                             ** dfAdult.meta.loc[dfAdult.meta.age >= 6,"selF"]).astype(np.int64)
################################################################
#    else: #base fecundity when no drugs
#        #all locations where age is less than 6
#        dfAdult.meta.loc[dfAdult.meta.age < 6, "fec"] = np.random.poisson(fecund,
#                  len(dfAdult.meta[dfAdult.meta.age < 6]))
#        #linear function defining decline in fecundity with age
#        m = float(0 - fecund) / (21 - 6)
#        b = 0 - m * 21
#        #assign fecundity value based on age function
#        positive_lambda = (dfAdult.meta.loc[dfAdult.meta.age >= 6, "age"].values * m) + b
#        positive_lambda[positive_lambda < 0] = 0
#        dfAdult.meta.loc[dfAdult.meta.age >= 6, "fec"] = np.random.poisson(positive_lambda).astype(np.int64)
#
#    #sex, recombination, mutation
#    dfAdult_mf = recombination_fx(locus, dfAdult, recombination_rate, basepairs)
#    # Positions is just the new positions
#    dfAdult_mf, new_positions = mutation_fx(locus, dfAdult_mf,
#         mutation_rate, recombination_rate, basepairs)
#    if selection: #dfAdult.sel will be updated here to same length as dfAdult_mf.pos
#        dfAdult_mf, dfAdult = selection_fx(dfAdult, dfAdult_mf, new_positions)
#
#    return(dfAdult_mf, dfAdult)
#
#def fecunditymda_sel2_fx(fecund,
#                    dfAdult,
#                    locus,
#                    mutation_rate,
#                    recombination_rate,
#                    basepairs,
#                    selection,
#                    densitydep_fec,
#                    mda_sterile,
#                    clear_count,
#                    mda_clear,
#                    dfHost):
#    '''function for reduced fecundity under mda option 2
#    option 2 is when the mutant are less fit then the wildtype when no mda
#    is being applied.
#    conditions: mda=True, selection=True, 2
#
#    Parameters
#    ----------
#    fecund: int
#         rate of the poisson distribution for offspring number
#    sterile_p: float
#         proportion of the adult parasite that are permanently sterilized from drugs
#    clear_time: int
#         number of months the drugs is effective or present in the host
#    clear_count: int
#         count of months since drug administered
#    dfAdult: dataframe
#         pandas dataframe containing adult parasites
#    dfMF: dataframe
#          pandas dataframe containing larval stage parasites
#
#    Returns
#    ------
#    dfAdult
#    dfMF
#    '''
#    if clear_count == 1: #permanent sterility
#         for index, row in dfHost[dfHost.MDA == 1].iterrows():
#             #randomly select dfAdults, change sex to "S" for sterile
#             try:
#                 mdarand = np.random.random(len(dfAdult.meta.hostidx == row.hostidx))
#                 mdasterile = dfHost.meta.ix[np.where(mdarand < mda_sterile
#                              ** dfAdult.meta[dfAdult.meta.hostidx == row.hostidx]["selF"])].index.values
#                 dfAdult.meta.ix[mdasterile, "sex"] = "S"
#             except ValueError:
#                 pass
#
#    if clear_count > 0 and clear_count <= mda_clear: #Drugs cause temporary sterility over clear_time
#         #linear function defining fecundity during drug clearance
#         mmda = float(fecund - 1) / (mda_clear - 1 )
#         bmda = 1 - mmda * 1
#         #new base fecundity under drugs
#         sterile_t = (mmda * clear_count + bmda)
#         assert sterile_t > 0
#         #assign value to dfAdult.fec
#         dfAdult.meta.loc[dfAdult.meta.age < 6, "fec"] = np.random.poisson(sterile_t **
#                    dfAdult.meta.loc[dfAdult.meta.age < 6, "selF"],
#                               len(dfAdult.meta[dfAdult.meta.age < 6]))
#
#        #linear function defining decline in fecundity with age
#         m = float(0 - sterile_t) / (21 - 6)
#         b = 0 - m * 21
#         #assign fecundity value based on age function
#         positive_lambda = (dfAdult.meta.loc[dfAdult.meta.age >= 6, "age"].values * m) + b
#         positive_lambda[positive_lambda < 0] = 0
#         dfAdult.meta.loc[dfAdult.meta.age >= 6, "fec"] = np.random.poisson((positive_lambda)
#                             ** dfAdult.meta.loc[dfAdult.meta.age >= 6,"selF"]).astype(np.int64)
#############################################
#    else: #base fecundity when no drugs
#         #assign value to dfAdult.fec
#         dfAdult.meta.loc[dfAdult.meta.age < 6, "fec"] = np.random.poisson(fecund **
#                    (1 - abs(1 - dfAdult.meta.loc[dfAdult.meta.age < 6, "selF"])))
#
#         #linear function defining decline in fecundity with age
#         m = float(0 - fecund) / (21 - 6)
#         b = 0 - m * 21
#         #assign fecundity value based on age function
#         positive_lambda = (dfAdult.meta.loc[dfAdult.meta.age >= 6, "age"].values * m) + b
#         positive_lambda[positive_lambda < 0] = 0
#         dfAdult.meta.loc[dfAdult.meta.age >= 6, "fec"] = np.random.poisson(positive_lambda
#                         ** (1 - abs(1 - dfAdult.meta.loc[dfAdult.meta.age >= 6,"selF"])),
#                            len(dfAdult.meta[dfAdult.meta.age >= 6]))
#############################################
#    #sex, recombination, mutation
#    dfAdult_mf = recombination_fx(locus, dfAdult, recombination_rate, basepairs)
#    # Positions is just the new positions
#    dfAdult_mf, new_positions = mutation_fx(locus, dfAdult_mf,
#         mutation_rate, recombination_rate, basepairs)
#    if selection: #dfAdult.sel will be updated here to same length as dfAdult_mf.pos
#        dfAdult_mf, dfAdult = selection_fx(dfAdult, dfAdult_mf, new_positions)
#
#    return(dfAdult_mf, dfAdult)
