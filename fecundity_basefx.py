# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 12:19:29 2017
@author: stsmall
"""
import numpy as np
import pandas
from scipy.stats import weibull_min

dfAdult = pd.DataFrame(
            {"village" : [],
             "host" : [],
             "age" : [],
             "R0net" : [],
             "hap" : [],
             "dip1" : [],
             "dip2" : []
             })
dfJuv = pd.DataFrame(
            {"village" : [],
             "host" : [],
             "age" : [],
             "R0net" : [],             
             "hap" : [],
             "dip1" : [],
             "dip2" : []
             })
dfMF = pd.DataFrame(
            {"village" : [],
             "host" : [],
             "age" : [],
             "R0net" : [],             
             "hap" : [],
             "dip1" : [],
             "dip2" : []
             })

###conditions: mda=False, selection=False

def fecundity_basefx(fecund):
    '''base fecundity function, simpliest scenario'''
    #conditions mda=False, selection=False      

    ageltsixIDX = which(dfAdult.age < 6)
    dfAdult[ageltsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(fecund))

    agegtsixIDX = which(dfAdult.age >= 6)        
    m = float(0 - fecund) / (21 - 6)
    b = 0 - m * 21    
    dfAdult[agegtsixIDX] = dfAdult.assign(fec=lambda dfAdult: 
        np.random.poisson(m * dfAdult.age[agegtsixIDX] + b))
    
    dfMF = dfMF.append(dfAdult * dfAdult.fecundity_base) #basically copy the adults 
      #values to dfMF the number times in fecundity column
        
if __name__ == '__main__':
    fecundity_basefx(fecund=20)


    