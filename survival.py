#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
from scipy.stats import weibull_min
from fecundity import fecunditybase_fx
from host_migration import hostmigration_fx
import random
   
def survivalbase_fx(month, surv_Juv, shapeMF, scaleMF, shapeAdult,
                    scaleAdult, dfMF, dfAdult, dfJuv, dfHost,
                    fecund, locus, mutation_rate, recombination_rate, 
                    basepairs, selection, hostmigrate):
   #(1, 0.866, 3.3, 10, 3.8, 8, dfMF, dfAdult, dfJuv, dfHost)  
                        
    '''base survival function
    Parameters
    ---------
    month: int
         time in months of the simulation
    surv_Juv: float
         survival prob of juvenille life stage
    shapeMF: float
         shape parameter for weibull distribution of MF
    scaleMF: int
         scale parameter for weibull distribution of MF
    shapeAdult: float
         shape parameter for weibull distribution of Adults
    scaleAdult: int
         scale parameter for weibull distribution of MF

    Returns
    -------
    dfMF
    dfAdult
    dfJuv
    dfHost
    
    '''
    #adult worms and hosts are only evaluated per year
    if month%12 == 0:
        #Adult survival is based on weibull cdf
        surv_adultrand = np.random.random(len(dfAdult))    
        surv_adultfxage = weibull_min.cdf(dfAdult.age, shapeAdult,loc=0,scale=scaleAdult)
        surviveAdult = np.where(surv_adultrand <= (1 - surv_adultfxage))
        dfAdult = dfAdult.iloc[surviveAdult]
        dfAdult.age = dfAdult.age + 1 #2 - 21
        
        ##host survival is from act table
        dfHost = dfHost[dfHost.age < dfHost.agedeath]
        #remove all worms with dead host.hostidx from all dataframes
        dfAdult = dfAdult.loc[dfAdult["hostidx"].isin(dfHost.hostidx)]
        dfJuv = dfJuv.loc[dfJuv["hostidx"].isin(dfHost.hostidx)]
        dfMF = dfMF.loc[dfMF["hostidx"].isin(dfHost.hostidx)]
        #add 1 year to all ages of hosts
        dfHost.age = dfHost.age + 1
        dfHost = hostmigration_fx(dfHost, hostmigrate)
        
    ##Juv is exponential 0.866; surv_Juv
    #dont include age 0 which just moved from transmission fx
    surv_juvrand = np.random.random(len(np.where(dfJuv.age > 0))[0])
    surviveJuv = np.where(surv_juvrand <= surv_Juv)
    dfJuv = dfJuv.iloc[surviveJuv]
    dfJuv.age = dfJuv.age + 1 # 1 - 13

    ##MF is weibull cdf
    surv_mfrand = np.random.random(len(dfMF))    
    surv_mffxage = weibull_min.cdf(dfMF.age,shapeMF,loc=0,scale=scaleMF)
    surviveMF = np.where(surv_mfrand <= (1 - surv_mffxage))
    dfMF = dfMF.iloc[surviveMF]
    dfMF.age = dfMF.age + 1 #2 - 12
    dfMF = dfMF[dfMF.age < 13] #hard cutoff at 12 months

    ##move Juv age 13 to adult age 1
    #dfJuv_new = pd.DataFrame({})
    dfJuv_new = dfJuv[dfJuv.age > 12]
    #reset age to adult   
    dfJuv_new.age = 1
    #increase R0net for next gen
    dfJuv_new.R0net = dfJuv_new.R0net + 1
    #append to adults
    dfAdult = dfAdult.append(dfJuv_new, ignore_index=True)
    #remove Juv age 13 from dfJuv
    dfJuv = dfJuv[dfJuv.age <= 12]
     
    ##call to fecundity fx to deepcopy adult to dfMF age 1
    #fecundity calls mutation/recombination
    dfAdult_mf = fecunditybase_fx(fecund, dfAdult, locus, mutation_rate, recombination_rate, basepairs)
    dfAdult_mf.age = 1
    dfAdult_mf.fec = 0
    dfAdult_mf.sex = [random.choice("MF") for i in range(len(dfAdult_mf))]
    dfMF = dfMF.append(dfAdult_mf, ignore_index=True)     

    return dfAdult, dfHost, dfJuv, dfMF if month%12 == 0 else dfJuv, dfMF
