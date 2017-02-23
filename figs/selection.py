#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import copy
def fitness_fx(locus,
               dfAdult_mf,
               dfAdult):
    ''' calculates mean fitness for each individual by summing fitness effects
    from dfSel for each position across all loci
    Parameters
    ---------
    locus : int
        number of loci
    dfAdult_mf : df
         data frame of adult worms containing genotype information
    dfSel : df
         data fram of fitness benefit for each allele
    Returns
    ------
    dfAdult_mf : df
         updated df for mf
    '''
    avg_over = len(dfAdult_mf.h2.keys())
    ninds = dfAdult_mf.meta.shape[0]
    fitF_ind = np.zeros(ninds)
    fitS_ind = np.zeros(ninds)

    for locus in dfAdult_mf.h2.keys():
        #these need to be the same length for the dot mult
#        assert dfAdult_mf.h1[locus].shape == dfAdult.sel[locus + "F"].shape
#        assert dfAdult_mf.h2[locus].shape == dfAdult.sel[locus + "F"].shape
#        assert dfAdult_mf.h1[locus].shape == dfAdult.sel[locus + "S"].shape
#        assert dfAdult_mf.h2[locus].shape == dfAdult.sel[locus + "S"].shape
        sum_selsites_S = np.dot(dfAdult_mf.h1[locus], dfAdult.sel[locus + "S"]) \
            + np.dot(dfAdult_mf.h2[locus], dfAdult.sel[locus + "S"])
        sum_selsites_F = np.dot(dfAdult_mf.h1[locus], dfAdult.sel[locus + "F"]) \
            + np.dot(dfAdult_mf.h2[locus], dfAdult.sel[locus + "F"])
####
        intsites_S = copy.copy(dfAdult.sel[locus + "S"])
        intsites_S[intsites_S > 0] = 1
        intsites_F = copy.copy(dfAdult.sel[locus + "F"])
        intsites_F[intsites_F > 0] = 1
        cds_sites_S = np.dot(dfAdult_mf.h1[locus], intsites_S) \
            + np.dot(dfAdult_mf.h2[locus], intsites_S)
        cds_sites_F = np.dot(dfAdult_mf.h1[locus], intsites_F) \
            + np.dot(dfAdult_mf.h2[locus], intsites_F)
####
        fitS_ind += (( (dfAdult.sel[locus + "St"] * 2) - cds_sites_S) + sum_selsites_S) / (dfAdult.sel[locus + "St"] * 2)
        fitF_ind += (( (dfAdult.sel[locus + "Ft"] * 2) - cds_sites_F) + sum_selsites_F) / (dfAdult.sel[locus + "Ft"] * 2)
####
    dfAdult_mf["fitS"] = fitS_ind / avg_over
    dfAdult_mf["fitF"] = fitF_ind / avg_over
    return(dfAdult_mf)

def selection_fx(dfAdult,
                 dfAdult_mf,
                 new_positions,
                 locus):
    '''recalculates DFE for new mutations and phenotype for new mf
    Parameters
    ---------
    dfSel : df
         updates this dataframe
    dfAdult_mf : df
         adds phenotype info
    dfMuts : df
         gives positions of new mutations for addition to dfSel
    locus : int
         number of loci
    Returns
    ------
    dfSel : df
         now updates with new mutations
    dfAdult_mf : df
         updated with phenotype
    '''
    for loc in dfAdult_mf.h2.keys(): #since this wont include 0
        selS = []
        selF = []
        for pos in new_positions[loc]: #this is the dict of positions
            if any([i <= pos <= j for i,j in dfAdult.coord[loc + "F"]]):
                 #shape = 4, mean = 1, scale = mean/shape
                 #here mean is mean_fitness, wildtype is assumed to be 1
                 selF.append(np.random.gamma(4, scale=0.25))
                 selS.append(0)
            elif any([i <= pos <= j for i,j in dfAdult.coord[loc + "S"]]):
                     selS.append(np.random.gamma(4, scale=0.25))
                     selF.append(0)
            else: #not in a cds
                selS.append(0)
                selF.append(0)
        dfAdult.sel[loc + "S"] = np.insert(dfAdult.sel[loc + "S"], new_positions[loc], selS)
        dfAdult.sel[loc + "F"] = np.insert(dfAdult.sel[loc + "F"], new_positions[loc], selF)
    dfAdult_mf = fitness_fx(locus, dfAdult_mf, dfAdult)
    return(dfAdult_mf, dfAdult)
