#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import random
import pandas as pd
def fitness_fx(locus,
               dfAdult_mf,
               dfSel):
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

    fitS = []
    fitF = []
#    print dfSel.head()
    for index, row in dfAdult_mf.iterrows():
        fitS_ind = []
        fitF_ind = []
#        print row["locus_" + str(1) + "_h1"]
#        print row["locus_" + str(1) + "_h2"]
        for loc in range(1,locus):
            fitS_ind.extend(dfSel.loc[dfSel["position"].isin(row.ix
                                        ["locus_" + str(loc) + "_h1"])]
                                        ['selS'][dfSel["locus"] == loc])
            fitF_ind.extend(dfSel.loc[dfSel["position"].isin(row.ix
                                        ["locus_" + str(loc) + "_h1"])]
                                        ['selF'][dfSel["locus"] == loc])
        fitS.append(round(np.mean(fitS_ind), 5))
        fitF.append(round(np.mean(fitF_ind), 5))
    dfAdult_mf["fitS"] = fitS
    dfAdult_mf["fitF"] = fitF
#    print(dfAdult_mf.head())
    return(dfAdult_mf)

def selection_fx(dfAdult_mf,
                 positions,
                 dfSel,
                 locus,
                 cds_coordinates):
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

    for loc in range(len(positions)):
        cds_positions = []
        muts_counter = []
        for start, end in cds_coordinates[loc]:
            cds_positions.extend([pos for pos in positions[loc] if pos >= start and pos <= end])
        muts_counter.append(cds_positions)
    newposlist = []
    for loc in range(len(positions)):
        locu = loc + 1
        for pos in muts_counter[loc]:
            if random.choice("SF") is "F":
                #shape = 4, mean = 1, scale = mean/shape
                #here mean is mean_fitness, wildtype is assumed to be 1
                selF = np.random.gamma(4, scale=0.25)
                selS = 1
            else:
                selS = np.random.gamma(4, scale=0.25)
                selF = 1
            newposlist.append([locu, pos, selS, selF])
    dfSel = pd.concat([dfSel, pd.DataFrame(newposlist, columns=dfSel.columns)],ignore_index=True)

    dfAdult_mf = fitness_fx(locus, dfAdult_mf, dfSel)
    return(dfAdult_mf, dfSel)
