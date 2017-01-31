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
from fecundity import fecundity_fx
from host_migration import hostmigration_fx
import random
   
def survivalmda_fx(month, surv_Juv, shapeMF, scaleMF, shapeAdult,
                    scaleAdult, dfMF, dfAdult, dfJuv, dfHost,
                    fecund, locus, mutation_rate, recombination_rate, 
                    basepairs, selection, hostmigrate, mda, mda_start, mda_num,
                    mda_freq, mda_coverage, mda_macro, mda_juvicide, mda_micro, 
                    mda_sterile, mda_clear):
     
     
    mda = [False, False]
    mda_start = [12, 12]
    mda_num = [6, 6] #how many mdas
    mda_freq = 12 #every 12 months
    mda_coverage = [0.8, 0.7] 
    mda_macro = 0.05
    mda_juvicide = 0.45
    mda_micro = 0.95
    mda_sterile = 0.35
    mda_clear = 6
     
     
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
    mda : boolean
         mda in village or not
    mda_start : list, int
         time in months that mda began in each village
    mda_num : list, int
         how many mdas to simulate
    mda_freq : int
         how often
    mda_coverage : float
         what percentage of the population get mda
    mda_macro: float
         percent of adults killed by drug
    mda_micro: float
         percent of MF killed by drug
    mda_juvcide: float
         percent of juveniile killed by drug; could be 0 or avg micro/macro
    mda_sterile : float
         percent of adult worms perm sterilized
    mda_clear : int
         clearance time
   
    Returns
    -------
    dfMF
    dfAdult
    dfJuv
    dfHost
    
    '''    
freq = 12
mda_start = 12
mda_clear = 6
mda_num = 6  
for month in range(120):
    if month < mda_start:
        print "no mda"
    elif month >= mda_start and month < (mda_start + freq * num + freq):    
         if month%freq == 0:
              clear_count = 1
              print "month %i, clear count is 1" %(month)
         elif (month - mda_clear)%freq == 0:
              clear_count = mda_clear
              print "month %i, clear count is %i" %(month, mda_clear)
         #else:
         #     clear_count = (month-mda_clear)%freq - mda_clear + 1
         #     print "month %i, clear count is %i" %(month, clear_count)
    elif month > (mda_start + freq * num):    
        print "mda over"
        clear_count = 0
        

    
    
    
    
    if clear_count == 1:
        ##MF
        mfkill = np.random.random(1,len(dfMF.age)) #array of random numbers
        if mfkill <= microcide: #what index in the array is less than microcide, 0.90
             dfMF.drop(dfMF[[mfkill]]) #remove rows
        ##Juv
        juvkill = np.random.random(1,len(dfJuv.age)) #array of random numbers
        if juvkill <= juvcide: #what index in the array is less than juvcide, 0.45
             dfJuv.drop(dfJuv[[juvkill]]) #remove rows        
        ##Adult
        adultkill = np.random.random(1,len(dfAdult.age)) #array of random numbers
        if adultkill <= macrocide: #what index in the array is less than macrocide, 0.05
             dfAdult.drop(dfAdult[[adultkill]]) #remove rows
    else:
         ##Juv is exponential 0.866   
         survjuv = np.random.random(len(dfJuv.age)) #array of random numbers
         killjuv = which(survjuv >= surv_Juv)    #compare random numbers to survival
         dfJuv.drop(dfJuv[[killjuv]])   #remove entire row from df if dies
         
         ##MF is weibull cdf
         survmf = np.random.random(len(dfMF.age)) #array of random numbers    
         surv_MF = apply(dfMF.age, 2, function(x) weibull_min.cdf(dfMF.age,shapeMF,
                         loc=0,scale=scaleMF)) #random number from weibull
         killmf = which(survmf <= surv_MF)     #compare random numbers
         dfMF.drop(dfMF[[killmf]]) #remove rows of MF
     
         #adult worms are only evaluated per year
         if month%12 == 0:
             #Adult is weibull cdf
             survadult = np.random.random(len(dfAdult.age)) #array of random numbers    
             surv_Adult = apply(dfAdult.age, 2, function(x) weibull_min.cdf(dfAdult.age,
                                shapeAdult,loc=0,scale=scaleAdult)) #weibull
             killadult = which(survadult <= surv_Adult)  #compare
             dfAdult.drop(dfAdult[[killadult]])        #remove row

    return dfMF, dfJuv, dfAdult

def survivalmda_sel1_fx(month=1, 
                   macrocide=0.05, 
                   microcide=0.90,
                   juvcide=0.45, 
                   clear_count=1, 
                   surv_Juv=0.866, 
                   shapeMF=3.3, 
                   scaleMF=12, 
                   shapeAdult=3.8, 
                   scaleAdult=8,
                   dfMF,
                   dfAdult,
                   dfJuv):
 
    '''base survival function
    Parameters
    ---------
    month: int
         time in months of the simulation
    macrocide: float
         percent of adults killed by drug
    microcide: float
         percent of MF killed by drug
    juvcide: float
         percent of juveniile killed by drug; could be 0 or avg micro/macro
    clear_count: int
         time since MDA, drug was administered
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
    
    '''
    
    if clear_count == 1:
        ##MF
        mfkill = np.random.random(1,len(dfMF.age)) #array of random numbers
        if mfkill <= microcide ** (dfMF.selS): #what index in the array is less than microcide, 0.90
             dfMF.drop(dfMF[[mfkill]]) #remove rows
        ##Juv
        juvkill = np.random.random(1,len(dfJuv.age)) #array of random numbers
        if juvkill <= juvcide ** (dfJuv.selS): #what index in the array is less than juvcide, 0.45
             dfJuv.drop(dfJuv[[juvkill]]) #remove rows        
        ##Adult
        adultkill = np.random.random(1,len(dfAdult.age)) #array of random numbers
        if adultkill <= macrocide ** (dfAdult.selS): #what index in the array is less than macrocide, 0.05
             dfAdult.drop(dfAdult[[adultkill]]) #remove rows        
    else:
         ##Juv is exponential 0.866   
         survjuv = np.random.random(len(dfJuv.age)) #array of random numbers
         killjuv = which(survjuv >= surv_Juv)    #compare random numbers to survival
         dfJuv.drop(dfJuv[[killjuv]])   #remove entire row from df if dies
         
         ##MF is weibull cdf
         survmf = np.random.random(len(dfMF.age)) #array of random numbers    
         surv_MF = apply(dfMF.age, 2, function(x) weibull_min.cdf(dfMF.age,shapeMF,
                         loc=0,scale=scaleMF)) #random number from weibull
         killmf = which(survmf <= surv_MF)     #compare random numbers
         dfMF.drop(dfMF[[killmf]]) #remove rows of MF
     
         #adult worms are only evaluated per year
         if month%12 == 0:
             #Adult is weibull cdf
             survadult = np.random.random(len(dfAdult.age)) #array of random numbers    
             surv_Adult = apply(dfAdult.age, 2, function(x) weibull_min.cdf(dfAdult.age,
                                shapeAdult,loc=0,scale=scaleAdult)) #weibull
             killadult = which(survadult <= surv_Adult)  #compare
             dfAdult.drop(dfAdult[[killadult]])        #remove row

    return dfMF, dfJuv, dfAdult
    
def survivalmda_sel2_fx(month=1, 
                   macrocide=0.05, 
                   microcide=0.90,
                   juvcide=0.45, 
                   clear_count=1, 
                   surv_Juv=0.866, 
                   shapeMF=3.3, 
                   scaleMF=12, 
                   shapeAdult=3.8, 
                   scaleAdult=8,
                   dfMF,
                   dfAdult,
                   dfJuv):
 
    '''base survival function
    Parameters
    ---------
    month: int
         time in months of the simulation
    macrocide: float
         percent of adults killed by drug
    microcide: float
         percent of MF killed by drug
    juvcide: float
         percent of juveniile killed by drug; could be 0 or avg micro/macro
    clear_count: int
         time since MDA, drug was administered
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
    
    '''
    
    if clear_count == 1:
        ##MF
        mfkill = np.random.random(1,len(dfMF.age)) #array of random numbers
        if mfkill <= microcide ** (dfMF.selS): #what index in the array is less than microcide, 0.90
             dfMF.drop(dfMF[[mfkill]]) #remove rows
        ##Juv
        juvkill = np.random.random(1,len(dfJuv.age)) #array of random numbers
        if juvkill <= juvcide ** (dfJuv.selS): #what index in the array is less than juvcide, 0.45
             dfJuv.drop(dfJuv[[juvkill]]) #remove rows        
        ##Adult
        adultkill = np.random.random(1,len(dfAdult.age)) #array of random numbers
        if adultkill <= macrocide ** (dfAdult.selS): #what index in the array is less than macrocide, 0.05
             dfAdult.drop(dfAdult[[adultkill]]) #remove rows        
    else: #(1-abs(avg_fitness-selS))
        
         ##Juv is exponential 0.866   
         avg_Juvfit = mean(dfJuv.selS) 
         survjuv = np.random.random(len(dfJuv.age)) #array of random numbers
         killjuv = which(survjuv >= (surv_Juv ** (1-abs(avg_Juvfit-selS))))    #compare random numbers to survival
         dfJuv.drop(dfJuv[[killjuv]])   #remove entire row from df if dies
         
         ##MF is weibull cdf
         avg_MFf = mean(dfMF.selS)
         survmf = np.random.random(len(dfMF.age)) #array of random numbers    
         surv_MF = apply(dfMF.age, 2, function(x) weibull_min.cdf(dfMF.age,shapeMF,
                         loc=0,scale=scaleMF)) #random number from weibull
         killmf = which(survmf <= (surv_MF ** (1-abs(avg_MFfit-selS))))      #compare random numbers
         dfMF.drop(dfMF[[killmf]]) #remove rows of MF
     
         #adult worms are only evaluated per year
         if month%12 == 0:
             avg_Adultfit = mean(dfAdult.selS)
             #Adult is weibull cdf
             survadult = np.random.random(len(dfAdult.age)) #array of random numbers    
             surv_Adult = apply(dfAdult.age, 2, function(x) weibull_min.cdf(dfAdult.age,
                                shapeAdult,loc=0,scale=scaleAdult)) #weibull
             killadult = which(survadult <= surv_Adult ** (1-abs(avg_Adultfit-selS)))  #compare
             dfAdult.drop(dfAdult[[killadult]])        #remove row

    return dfMF, dfJuv, dfAdult         
        
if __name__ == '__main__':
#    survival_basefx(surv_Juv=0.866, shapeMF=3.3, scaleMF=12, shapeAdult=3.8, scaleAdult=8)
#    survival_mdafx(month=1, macrocide=0.05, microcide=0.90, juvcide-0.45, clear_count=1, 
#                   surv_Juv=0.866, shapeMF=3.3, scaleMF=12, shapeAdult=3.8, scaleAdult=8)
#    survival_sel1fx(month=1, macrocide=0.05, microcide=0.90, juvcide-0.45, clear_count=1, 
#                   surv_Juv=0.866, shapeMF=3.3, scaleMF=12, shapeAdult=3.8, scaleAdult=8)
#    survival_sel2fx(month=1, macrocide=0.05, microcide=0.90, juvcide-0.45, clear_count=1, 
#                   surv_Juv=0.866, shapeMF=3.3, scaleMF=12, shapeAdult=3.8, scaleAdult=8)
