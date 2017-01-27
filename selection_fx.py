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

def sel_fx(locus, positions):
     '''initializes the distribution of fitness effects for each mutation
     
     Parameters
     ---------
     locus : int
          number of loci
     positions : list
          list from ms_outcall, the the line "positions" in the scrm/ms output
     Returns
     ------
     dfSel
     '''
     selF =[]
     selS =[]
     positions = [item for sublist in positions for item in sublist]
     #below functs assume the list is nested 
     if isinstance(positions[0], list):
          numpos = [len(i) for i in positions]
          for loc in range(locus):
               for pos in positions[loc]:
                    if random.choice("SF") is "F":
                         #shape = 4, mean = 1, scale = mean/shape
                         #here mean is mean_fitness, wildtype is assumed to be 1
                         selF.append(np.random.gamma(4, scale=0.25))
                         selS.append(1)
                    else:
                         selS.append(np.random.gamma(4, scale=0.25))
                         selF.append(1)
         
          dfSel = pd.DataFrame({
                                'locus' : np.repeat(range(1, locus), numpos),
                                'position' : [item for sub in positions for item in sub],
                                'selF' : selF,
                                'selS' : selS,
                                'freqInt' : np.zeros(sum(numpos))})     
          dfSel = dfSel.loc[:, ['locus', 'position', 'selF',
                 'selS', 'freqInit']]
     else: #list is not nested
          numpos = len(positions)        
          for pos in positions:
               if random.choice("SF") is "F":
                    #shape = 4, mean = 1, scale = mean/shape
                    #here mean is mean_fitness, wildtype is assumed to be 1
                    selF.append(np.random.gamma(4, scale=0.25))
                    selS.append(1)
               else:
                    selS.append(np.random.gamma(4, scale=0.25))
                    selF.append(1)
         
          dfSel = pd.DataFrame({
                                'locus' : np.repeat(range(1, locus), numpos),
                                'position' : [item for item in positions],
                                'selF' : selF,
                                'selS' : selS,
                                'freqInt' : np.zeros(numpos)})     
          dfSel = dfSel.loc[:, ['locus', 'position', 'selF',
                 'selS', 'freqInit']]  
     return dfSel
     
def fitness_fx(locus, dfAdult, dfSel):
     ''' calculates mean fitness for each individual by summing fitness effects
     from dfSel for each position across all loci
     
     Parameters
     ---------
     dfAdult : df
          data frame of adult worms containing genotype information
     dfSel : df
          data fram of fitness benefit for each allele
     
     Returns
     ------
     fitF : array
          array filling selF column, influencing fecundity
     fitS : array      
          array filling selS column, influencing survival
     freqinit : array
          initial freq of mutations for dfSel
     '''
     allele = []
    
     ##freq of positions in dfSel
     for loc in range(1, locus):
          for row in range(len(dfAdult)):
               #flatten list and extend
               allele.extend(dfAdult.iloc[row]["locus_" + str(loc) + "_h1"])
               allele.extend(dfAdult.iloc[row]["locus_" + str(loc) + "_h2"])
     freq = np.unique(allele, return_counts=True)[1] / (len(dfAdult) * 2.0) 
     
     ##fitness of individual in dfAdult from values in dfSel               
     fitS_ind = []
     fitF_ind = []
     fitS = []
     fitF = []
     for row in range(len(dfAdult)):
          for loc in range(1, locus):                         
               fitS_ind.extend(dfSel.loc[dfSel["position"].isin(dfAdult.iloc[row]
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selS'][dfSel["locus"] == loc])
               fitF_ind.extend(dfSel.loc[dfSel["position"].isin(dfAdult.iloc[row]
                                         ["locus_" + str(loc) + "_h1"])]
                                         ['selF'][dfSel["locus"] == loc])
          fitS.append(round(np.mean(fitS_ind), 5))
          fitF.append(round(np.mean(fitF_ind), 5))  
    
     return fitS, fitF, freq 