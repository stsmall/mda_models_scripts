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

def allele_freq_fx(dfAdult, dfSel):
    ''' Calculates allele freq for mutations in dfSel by village
    Parameters
    ---------
    dfAdult : df
        data frame of adult worms
    dfSel : df
        data frame of mutation positions at each locus
    Returns
    -------
    dfFreq : list, float
        freq of mutations in dfSel for each village
    '''
    #list of things to count from dfSel per locus, position
    positions = dfSel.groupby("locus").position.apply(lambda x: x.values.tolist())
    #counts of each mutation in dfAdult by village for each locus
    for loc in range(1,len(positions)):
        h1 = dfAdult.groupby("village")["locus_" + str(loc) + "_h1"].apply(lambda x: np.unique([i for l in 
                       x.values.tolist() for i in l], return_counts=True))
        h2 = dfAdult.groupby("village")["locus_" + str(loc) + "_h2"].apply(lambda x: np.unique([i for l in 
                       x.values.tolist() for i in l], return_counts=True))
    #h1[0][0] # positions in village 1 hap1
    #h1[0][1] # counts of positions in village 1 hap1
    #h1[1][0] # positions in village 2 hap1
    #h1[1][1] # counts of positions in village 2 hap1
    dfFreq = pd.DataFrame({"village" : [],
                           "locus" : np.repeat,
                           "position" : [positions],
                           "allele_freq" : []             
                           })
    