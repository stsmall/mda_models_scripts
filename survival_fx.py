# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 12:19:29 2017
@author: stsmall
"""
import numpy as np
import pandas
from scipy.stats import weibull_min

def survival_host(dfHost, agebound=[5,80]):
    '''Host can die, the death of a host signifies the entire loss of a worm 
    population. Host ages are choosen using actuary tables, host deaths are 
    probabilities based on cumulative survival rates. Deaths are evaluated per year
    This function SHOULD be referenced from the parasite survival functions inside 
    the condition if 'month%12 is 0:'
    
    Parameters
    ---------
    dfHost: df
        dataframe to be updated if host dies
    agebound: int,list
        list of youngest and oldest ages allowed
    fertHost: int
        estimate of fertility for new age
    
    Returns
    -------
    dfHost
    '''
    young, old = 5, 80 #age bounds
    survCum = dfHost.apply(dfHost.age, 2, function(x) weibull_min.cdf(young,old))
    survHost = np.random.random(len(dfHost.age)) #array of random numbers
    killHost = which(survHost <= survCum)    #compare random numbers to survival
    dfHost.drop(dfHost[[killHost]])   #remove entire row from df if dies    
    
    ###fertHost=5
    #this should ideally go with new infection since that is when a new host is selected
    #newHost = np.random.geometric(mean, young, old) #age of new host
    
def survival_basefx(month=1, 
                    surv_Juv=0.866, 
                    shapeMF=3.3, 
                    scaleMF=10, 
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
        
        
def survival_mdafx(month=1, 
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

def survival_sel1fx(month=1, 
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
    
def survival_sel2fx(month=1, 
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
