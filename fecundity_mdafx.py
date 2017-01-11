# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 17:40:02 2017

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
             "fec" : [],      #main change from base fxs
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
### conditions mda=T, selection=F

def fecundity_mdafx(fecund=20, 
                    sterile_p=0.35, 
                    clear_time=6, 
                    clear_count=1):
    '''function for reduced fecundity under mda'''

    #permanent sterility from drugs
    if clear_count == 1: #immediately folowing treatment
        sterile = np.random.randint(1,len(dfAdult.age),
                                    round(sterile_p * len(dfAdult.age)))
        dfAdult.fec[[sterile]] = 0  #permanently sterlized adults
    
    #basefecundity, same as fecundity_basefx
    ageltsixIDX = which(dfAdult.age < 6 && dfAdult.fec != 0)
    dfAdult[ageltsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(fecund))

    agegtsixIDX = which(dfAdult.age >= 6  && dfAdult.fec != 0)        
    m = float(0 - fecund) / (21 - 6)
    b = 0 - m * 21    
    dfAdult[agegtsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(m * dfAdult.age[agegtsixIDX] + b))

    #temp sterility from drugs
    if clear_count < 6:        
        m = float(0.9 - 1) / (1 - clear_time)
        b = 0.9 - m * 1
        sterile_t = (m * clear_count + b)  
    else:
        sterile_t = 1

    #penalize fecundity for temp sterility        
    dfAdult.fec = round(dfAdult.fec * sterile_t)  
     
    #produce MF by copying adult row to dfMF 
    dfMF = dfMF.append(dfAdult * dfAdult.fec) #basically copy the adults 
      #values to dfMF the number times in fecundity column
    
    return dfAdult, dfMF

if __name__ == '__main__':
    fecundity_mdafx(fecund=20, sterile_p=0.35, clear_time=6, clear_count=1)
    