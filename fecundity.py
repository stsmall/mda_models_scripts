#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
from recombination import recombination_fx
from mutation import mutation_fx
from selection import selection_fx

def fecundity_basefx(fecund, dfAdult, locus, mutation_rate, recombination_rate, basepairs, selection, dfSel):
    #(20, dfAdult, 2, [7.6E-8, 2.9E-9], [0, 2.9E-9], [13000, 200000])
     '''base fecundity function, simpliest scenario
    conditions: mda=False, selection=F
    
    Parameters
    ---------
    fecund: int
         rate of the poisson distribution for offspring number
    dfAdult: pandas dataframe
          dataframe containing adult parasites
    dfMF: dataframe
          pandas dataframe containing larval stage parasites 
          
    Returns
    ------
    dfAdult_mf : df
         deep copy of adult genotypes
    
    TO DO: track the rate of decline for fecundity, should be linear
     
     '''
     #all locations where age is less than 6
     dfAdult.loc[dfAdult.age < 6, "fec"] = np.random.poisson(fecund, 
               len(dfAdult[dfAdult.age < 6]))         
     #linear function defining decline in fecundity with age
     m = float(0 - fecund) / (21 - 6)
     b = 0 - m * 21
     #assign fecundity value based on age function     
     dfAdult.loc[dfAdult.age >= 6, "fec"] = np.random.poisson(m 
               * dfAdult.loc[dfAdult.age >= 6,"age"] + b)

     #sex, recombination, mutation
     dfAdult_mf = recombination_fx(locus, dfAdult, recombination_rate, basepairs)
     dfAdult_mf, dfMuts = mutation_fx(locus, dfAdult_mf, mutation_rate, recombination_rate, basepairs)
    
     if selection:
         dfAdult_mf, dfSel = selection_fx(dfAdult_mf, dfMuts, dfSel, locus)
    
     return dfAdult_mf, dfSel if selection is True else dfAdult_mf
