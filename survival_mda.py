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
from fecundity_mda import fecunditymda_fx
from fecundity_mda import fecunditymda_sel1_fx
from fecundity_mda import fecunditymda_sel2_fx
from host_migration import hostmigration_fx
from hostmda import hostmda_fx
import random
   
def survivalmda_fx(month, villages, surv_Juv, shapeMF, scaleMF, shapeAdult,
                    scaleAdult, dfMF, dfAdult, dfJuv, dfHost,
                    fecund, locus, mutation_rate, recombination_rate, 
                    basepairs, selection, hostmigrate, mda, mda_start, mda_num,
                    mda_freq, mda_coverage, mda_macro, mda_juvicide, mda_micro, 
                    mda_sterile, mda_clear):
                            
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
#    mda = [False, False]
#    mda_start = [12, 12]
#    mda_num = [6, 6] #how many mdas
#    mda_freq = 12 #every 12 months
#    mda_coverage = [0.8, 0.7] 
#    mda_macro = 0.05
#    mda_juvicide = 0.45
#    mda_micro = 0.95
#    mda_sterile = 0.35
#    mda_clear = 6
     
    ##calculates time since mda as clear_count    
    if month < mda_start:
            clear_count = 0
    elif month >= mda_start and month < (mda_start + mda_freq * mda_num + mda_freq):    
        if month%mda_freq == 0:
             clear_count = 1
        elif (month - mda_clear)%mda_freq == 0:
             clear_count = mda_clear + 1
        else:
             if ((month-mda_freq)%mda_freq + 1) <= mda_clear:
                  clear_count = ((month-mda_freq)%mda_freq + 1)
             else:
                  clear_count = 0
    elif month > (mda_start + mda_freq * mda_num):    
        clear_count = 0
    
    ##assign MDA to hosts
    if clear_count == 1:
        ##add mda to dfHost per coverage    
        dfHost = hostmda_fx(villages, dfHost, mda_coverage)
        ##apply MDA effect to host populations of worms    
        for index, row in dfHost[dfHost.MDA == 1].iterrows():   
             #kill MF    
             dfMF.drop(dfMF[dfMF.hostidx == row.hostidx].sample(frac = mda_micro).index, inplace = True)    
             #kill Juv
             dfJuv.drop(dfJuv[dfJuv.hostidx == row.hostidx].sample(frac = mda_juvicide).index, inplace = True)
             #kill Adults
             dfAdult.drop(dfAdult[dfAdult.hostidx == row.hostidx].sample(frac = mda_macro).index, inplace = True)
    
    ##normal surival funxtions      
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
    dfAdult_mf = fecunditymda_fx(fecund, dfAdult, locus, mutation_rate, recombination_rate, 
                                 basepairs, clear_count, mda_sterile, mda_clear, dfHost, villages)
    dfAdult_mf.age = 1
    dfAdult_mf.fec = 0
    dfAdult_mf.sex = [random.choice("MF") for i in range(len(dfAdult_mf))]
    dfMF = dfMF.append(dfAdult_mf, ignore_index=True)     

    return dfAdult, dfHost, dfJuv, dfMF if month%12 == 0 else dfJuv, dfMF

def survivalmda_sel1_fx(month, villages, surv_Juv, shapeMF, scaleMF, shapeAdult,
                    scaleAdult, dfMF, dfAdult, dfJuv, dfHost,
                    fecund, locus, mutation_rate, recombination_rate, 
                    basepairs, selection, hostmigrate, mda, mda_start, mda_num,
                    mda_freq, mda_coverage, mda_macro, mda_juvicide, mda_micro, 
                    mda_sterile, mda_clear):

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
    
    ##calculates time since mda as clear_count    
    if month < mda_start:
            clear_count = 0
    elif month >= mda_start and month < (mda_start + mda_freq * mda_num + mda_freq):    
        if month%mda_freq == 0:
             clear_count = 1
        elif (month - mda_clear)%mda_freq == 0:
             clear_count = mda_clear + 1
        else:
             if ((month-mda_freq)%mda_freq + 1) <= mda_clear:
                  clear_count = ((month-mda_freq)%mda_freq + 1)
             else:
                  clear_count = 0
    elif month > (mda_start + mda_freq * mda_num):    
        clear_count = 0
    
    ##assign MDA to hosts
    if clear_count == 1:
        ##add mda to dfHost per coverage    
        dfHost = hostmda_fx(villages, dfHost, mda_coverage)
         ##apply MDA effect to host populations of worms    
        for index, row in dfHost[dfHost.MDA == 1].iterrows():   
             mdarand = np.random.random(len(dfMF.hostidx == row.hostidx))
             mdakill = np.where(mdarand < mda_micro ** dfMF[dfMF.hostidx == row.hostidx]["selS"])     
             dfMF.drop(dfMF.iloc[mdakill],inplace=True)

             mdarand = np.random.random(len(dfJuv.hostidx == row.hostidx))
             mdakill = np.where(mdarand < mda_juvicide ** dfJuv[dfJuv.hostidx == row.hostidx]["selS"])     
             dfJuv.drop(dfJuv.iloc[mdakill],inplace=True)

             mdarand = np.random.random(len(dfAdult.hostidx == row.hostidx))
             mdakill = np.where(mdarand < mda_macro ** dfAdult[dfAdult.hostidx == row.hostidx]["selS"])     
             dfAdult.drop(dfAdult.iloc[mdakill],inplace=True)             
    ##normal surival funxtions      
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
    dfAdult_mf = fecunditymda_sel1_fx(fecund, dfAdult, locus, mutation_rate, recombination_rate, 
                                 basepairs, clear_count, mda_sterile, mda_clear, dfHost, villages)
    dfAdult_mf.age = 1
    dfAdult_mf.fec = 0
    dfAdult_mf.sex = [random.choice("MF") for i in range(len(dfAdult_mf))]
    dfMF = dfMF.append(dfAdult_mf, ignore_index=True)     

    return dfAdult, dfHost, dfJuv, dfMF if month%12 == 0 else dfJuv, dfMF
    
def survivalmda_sel2_fx(month, villages, surv_Juv, shapeMF, scaleMF, shapeAdult,
                    scaleAdult, dfMF, dfAdult, dfJuv, dfHost,
                    fecund, locus, mutation_rate, recombination_rate, 
                    basepairs, selection, hostmigrate, mda, mda_start, mda_num,
                    mda_freq, mda_coverage, mda_macro, mda_juvicide, mda_micro, 
                    mda_sterile, mda_clear):
 
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
    
    ##calculates time since mda as clear_count    
    if month < mda_start:
            clear_count = 0
    elif month >= mda_start and month < (mda_start + mda_freq * mda_num + mda_freq):    
        if month%mda_freq == 0:
             clear_count = 1
        elif (month - mda_clear)%mda_freq == 0:
             clear_count = mda_clear + 1
        else:
             if ((month-mda_freq)%mda_freq + 1) <= mda_clear:
                  clear_count = ((month-mda_freq)%mda_freq + 1)
             else:
                  clear_count = 0
    elif month > (mda_start + mda_freq * mda_num):    
        clear_count = 0
    
    ##assign MDA to hosts
    if clear_count == 1:
        ##add mda to dfHost per coverage    
        dfHost = hostmda_fx(villages, dfHost, mda_coverage)
         ##apply MDA effect to host populations of worms    
        for index, row in dfHost[dfHost.MDA == 1].iterrows():   
             mdarand = np.random.random(len(dfMF.hostidx == row.hostidx))
             mdakill = np.where(mdarand < mda_micro ** dfMF[dfMF.hostidx == row.hostidx]["selS"])     
             dfMF.drop(dfMF.iloc[mdakill],inplace=True)

             mdarand = np.random.random(len(dfJuv.hostidx == row.hostidx))
             mdakill = np.where(mdarand < mda_juvicide ** dfJuv[dfJuv.hostidx == row.hostidx]["selS"])     
             dfJuv.drop(dfJuv.iloc[mdakill],inplace=True)

             mdarand = np.random.random(len(dfAdult.hostidx == row.hostidx))
             mdakill = np.where(mdarand < mda_macro ** dfAdult[dfAdult.hostidx == row.hostidx]["selS"])     
             dfAdult.drop(dfAdult.iloc[mdakill],inplace=True)             
    
    ##normal surival funxtions      
    #adult worms and hosts are only evaluated per year
    if month%12 == 0:
        #Adult survival is based on weibull cdf
        surv_adultrand = np.random.random(len(dfAdult))    
        surv_adultfxage = weibull_min.cdf(dfAdult.age, shapeAdult,loc=0,scale=scaleAdult)        
        surv_selfx = surv_adultfxage ** (1-abs(1-dfAdult.selS))        
        surviveAdult = np.where(surv_adultrand <= (1 - surv_selfx))
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
    surv_juvarray = np.repeat(surv_Juv, len(surv_juvrand))
    surv_selfx = surv_juvarray ** (1-abs(1-dfJuv.selS))
    surviveJuv = np.where(surv_juvrand <= surv_selfx)
    dfJuv = dfJuv.iloc[surviveJuv]
    dfJuv.age = dfJuv.age + 1 # 1 - 13

    ##MF is weibull cdf
    surv_mfrand = np.random.random(len(dfMF))    
    surv_mffxage = weibull_min.cdf(dfMF.age,shapeMF,loc=0,scale=scaleMF)
    surv_selfx = surv_mffxage ** (1-abs(1-dfMF.selS))
    surviveMF = np.where(surv_mfrand <= (1 - surv_selfx))
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
    dfAdult_mf = fecunditymda_sel2_fx(fecund, dfAdult, locus, mutation_rate, recombination_rate, 
                                 basepairs, clear_count, mda_sterile, mda_clear, dfHost, villages)
    dfAdult_mf.fec = 0
    dfAdult_mf.sex = [random.choice("MF") for i in range(len(dfAdult_mf))]
    dfMF = dfMF.append(dfAdult_mf, ignore_index=True)     

    return dfAdult, dfHost, dfJuv, dfMF if month%12 == 0 else dfJuv, dfMF