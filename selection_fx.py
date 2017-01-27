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
import random

def selection_fx(dfAdult_mf, dfMuts, dfSel, locus):
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
     for loc in range(1, locus):
          for index, row in dfMuts[dfMuts.locus == loc].iterrows():
               if row["position"].isin[dfSel[dfSel.locus == loc]]:
                    pass
               else:
                    if random.choice("SF") is "F":
                         #shape = 4, mean = 1, scale = mean/shape
                         #here mean is mean_fitness, wildtype is assumed to be 1
                         row.selF.set_value(np.random.gamma(4, scale=0.25))
                         row.selS.set_value(1)
                    else:
                         row.selS.set_value(np.random.gamma(4, scale=0.25))
                         row.selF.set_value(1)        
     dfSel = dfSel.append(dfMuts)
     dfAdult_mf = fitness_fx(locus, dfAdult_mf, dfSel)               
     return dfAdult_mf, dfSel
     
def fitness_fx(locus, dfAdult_mf, dfSel):
     ''' calculates mean fitness for each individual by summing fitness effects
     from dfSel for each position across all loci
     
     Parameters
     ---------
     dfAdult_mf : df
          data frame of adult worms containing genotype information
     dfSel : df
          data fram of fitness benefit for each allele
     
     Returns
     ------
     dfAdult_mf : df
          updated df for mf
     '''           
     for index, row in dfAdult_mf.iterrows():
          fitS_ind = []
          fitF_ind = []
          for loc in range(1,locus):   
               fitS_ind.extend(dfSel.loc[dfSel["position"].isin(row
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selS'][dfSel["locus"] == loc])
               fitF_ind.extend(dfSel.loc[dfSel["position"].isin(row
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selF'][dfSel["locus"] == loc])
          row.fitS.set_value(round(np.mean(fitS_ind), 5))
          row.fitF.set_value(round(np.mean(fitF_ind), 5))  
    
     return dfAdult_mf