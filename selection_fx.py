# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

def phenotype_fx(dfSel, df):
     '''sums fitness of phenotype by net value of mutation fitness
     Parameters
     --------
     dfSel: df
          dataframe of positions with phenotype contribution
     df: df
          dataframe of MF, since this is the only time new mutations happen
          once this is calculated it just copies to dfJuv and dfAdult
     
     Returns
     ------
     df, the dataframe input with the selF and selS columns totalled 
     '''
     
     return df
     
def DFE_fx(positions, seloption):
     '''initializes the distribution of fitness effects for each mutation
     
     Parameters
     ---------
     positions: list
          list from ms_outcall, the the line "positions" in the scrm/ms output
     seloption: int
          0 is no selection
          1 is no historic seleciton
          2 is historic selection
     Returns
     ------
     dfSel
     '''
     #draw from distribution of fitness effect for each mutation
     dfSel = pd.DataFrame({
                           'locus' : int,
                           'position' : int,
                           'selF' : float,
                           'selS' : float,
                           'freqInt' : float})
     
     return dfSel
     