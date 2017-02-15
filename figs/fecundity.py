#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np

from .recombination import recombination_fx
from .mutation import mutation_fx
from .selection import selection_fx

from IPython import embed


def fecunditybase_fx(fecund,
                     dfAdult,
                     locus,
                     mutation_rate,
                     recombination_rate,
                     basepairs,
                     selection,
                     dfSel,
                     cds_coordinates,
                     densitydep_fec):
    '''Base fecundity function, simpliest scenario

    Parameters
    ---------
    fecund: int
        rate of the poisson distribution for offspring number
    dfAdult: figs.Worms 
        dataframe containing adult parasites
    locus : int
    mutation_rate : list, float
    recombination_rate : list, float
    basepairs : list, int
    selection : boolean
    dfSel : df
    cds_coordinates : [list], int
    
    Returns
    ------
    dfAdult_mf : df
         deep copy of adult genotypes
    dfSel : df

    '''
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
    dfAdult_mf, positions = mutation_fx(locus, dfAdult_mf, 
         mutation_rate, recombination_rate, basepairs)

    if selection:
     dfAdult_mf, dfSel = selection_fx(dfAdult_mf, positions, dfSel, locus, cds_coordinates)

    return(dfAdult_mf, dfSel)
