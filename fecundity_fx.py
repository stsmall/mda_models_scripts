#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import pandas
from scipy.stats import weibull_min

def fecundity_basefx(fecund, dfAdult, locus, mutation_rate, recombination_rate, basepairs):
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
         dfAdult_mf, dfSel = DFE_fx(dfAdult_mf, dfMuts, dfSel)
    
     return dfAdult_mf, dfSel if selection is True else dfAdult_mf

####above code works

def fecundity_mdafx(fecund=20, 
                    sterile_p=0.35, 
                    clear_time=6, 
                    clear_count=1, 
                    dfAdult, 
                    dfMF):
    '''function for reduced fecundity under mda
    conditions: mda=True, selection=False
    
    Parameters
    ----------
    fecund: int
         rate of the poisson distribution for offspring number
    sterile_p: float
         proportion of the adult parasite that are permanently sterilized from drugs
    clear_time: int
         number of months the drugs is effective or present in the host
    clear_count: int
         count of months since drug administered
    dfAdult: dataframe
         pandas dataframe containing adult parasites
    dfMF: dataframe
          pandas dataframe containing larval stage parasites 
         
    Returns
    ------
    dfAdult
    dfMF
    '''
    if clear_count == 1: #permanent sterility
        sterile = np.random.random(1,len(dfAdult.age))
        if sterile <= sterile_p:
             #permanently sterilized adults
             dfAdult.fec[[sterile]] = 0
             
    elif clear_count < clear_time: #Drugs cause temporary sterility over clear_time
        #linear function defining fecundity during drug clearance
        m = float(fecund - 1) / (clear_time - 1 )
        b = 1 - m * 1
        #new base fecundity under drugs
        sterile_t = (m * clear_count + b)
        #assign value to dfAdult.fec
        ageltsixIDX = which(dfAdult.age < 6 && dfAdult.fec != 0)
        dfAdult[ageltsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(sterile_t))
        #assign value to dfAdult.fec with older worms
        agegtsixIDX = which(dfAdult.age >= 6  && dfAdult.fec != 0)
        #linear function of reproductive decline         
        mage = float(0 - fecund) / (21 - 6)
        bage = 0 - m * 21
        #sterile_t = ((float((mage * dfAdult.age[agegtsix] + bage)-1)/(clear_time-1)) 
         #    * clear_count + (1 - (float((mage * dfAdult.age[agegtsix] + bage)-1) 
          #        / (clear_time-1))))
        dfAdult[agegtsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(
             ((float((mage * dfAdult.age[agegtsix] + bage)-1)/(clear_time-1)) 
             * clear_count + (1 - (float((mage * dfAdult.age[agegtsix] + bage)-1) 
                  / (clear_time-1))))))
        
    else: #base fecundity when no drugs   
        #find all indexes where age is less than 6 and adult is not premanently sterilized
        ageltsixIDX = which(dfAdult.age < 6 && dfAdult.fec != 0)
        #assign value to dfAdult.fec
        dfAdult[ageltsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(fecund))
        #find all indexes (rows) where age is greater than equal to 6 and not sterilized
        agegtsixIDX = which(dfAdult.age >= 6  && dfAdult.fec != 0)
        #linear function of reproductive decline         
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign value to dfAdult.fec from linear function     
        dfAdult[agegtsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson((m 
              * dfAdult.age[agegtsixIDX] + b))
    
    #actual reproduction just copies adult info from dfAdult to dfMF
    #number of times indicated in the dfAdult.fec column
    dfMF = dfMF.append(rep(dfAdult,dfAdult.fec))

    return dfAdult, dfMF 
    
def fecundity_mdasel1fx(fecund=20, 
                    sterile_p=0.35, 
                    clear_time=6, 
                    clear_count=1, 
                    dfAdult, 
                    dfMF):
     
    '''function for reduced fecundity under mda option 1
    option 1 simplifies that when no MDA or selective event all phenotypes 
    are essetially wildtype, so fitness is not evaluated
    conditions: mda=True, selection=True, 1
    
    Parameters
    ----------
    fecund: int
         rate of the poisson distribution for offspring number
    sterile_p: float
         proportion of the adult parasite that are permanently sterilized from drugs
    clear_time: int
         number of months the drugs is effective or present in the host
    clear_count: int
         count of months since drug administered
    dfAdult: dataframe
         pandas dataframe containing adult parasites
    dfMF: dataframe
          pandas dataframe containing larval stage parasites     
         
    Returns
    ------
    dfAdult
    dfMF
    '''
    if clear_count == 1: #permanent sterility
        sterile = np.random.random(1,len(dfAdult.age))
        #here dfAdult.selF is the sum of allele effect on phenotype
        # 1 is wildtype, <1 is less fit than wt, >1 is more fit than wt
        if sterile <= sterile_p ** (dfAdult.selF):
             #permanently sterilized adults
             dfAdult.fec[[sterile]] = 0
             
    elif clear_count < clear_time: #Drugs cause temporary sterility over clear_time
        #linear function defining fecundity during drug clearance
        m = float(fecund - 1) / (clear_time - 1 )
        b = 1 - m * 1
        #new base fecundity under drugs
        sterile_t = (m * clear_count + b)
        #assign value to dfAdult.fec
        ageltsixIDX = which(dfAdult.age < 6 && dfAdult.fec != 0)
        dfAdult[ageltsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson((sterile_t ** dfAdult.selF)))

        #assign value to dfAdult.fec with older worms
        agegtsixIDX = which(dfAdult.age >= 6  && dfAdult.fec != 0)
        #linear function of reproductive decline         
        mage = float(0 - fecund) / (21 - 6)
        bage = 0 - m * 21
        #sterile_t = ((float((mage * dfAdult.age[agegtsix] + bage)-1)/(clear_time-1)) 
         #    * clear_count + (1 - (float((mage * dfAdult.age[agegtsix] + bage)-1) 
          #        / (clear_time-1))))
        dfAdult[agegtsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(
             (((float((mage * dfAdult.age[agegtsix] + bage)-1)/(clear_time-1)) 
             * clear_count + (1 - (float((mage * dfAdult.age[agegtsix] + bage)-1) 
                  / (clear_time-1)))) ** dfAdult.selF)))
    else: #base fecundity when no drugs   
        #find all indexes where age is less than 6 and adult is not premanently sterilized
        ageltsixIDX = which(dfAdult.age < 6 && dfAdult.fec != 0)
        #assign value to dfAdult.fec
        dfAdult[ageltsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(fecund))
        #find all indexes (rows) where age is greater than equal to 6 and not sterilized
        agegtsixIDX = which(dfAdult.age >= 6  && dfAdult.fec != 0)
        #linear function of reproductive decline         
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign value to dfAdult.fec from linear function     
        dfAdult[agegtsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson((m 
              * dfAdult.age[agegtsixIDX] + b)))
    
    #actual reproduction just copies adult info from dfAdult to dfMF
    #number of times indicated in the dfAdult.fec column
    dfMF = dfMF.append(rep(dfAdult,dfAdult.fec))

    return dfAdult, dfMF 
    
def fecundity_mdasel2fx(fecund=20, 
                    sterile_p=0.35, 
                    clear_time=6, 
                    clear_count=1, 
                    dfAdult, 
                    dfMF):
     
    '''function for reduced fecundity under mda option 2
    option 2 is when the mutant are less fit then the wildtype when no mda
    is being applied.
    conditions: mda=True, selection=True, 2
    
    Parameters
    ----------
    fecund: int
         rate of the poisson distribution for offspring number
    sterile_p: float
         proportion of the adult parasite that are permanently sterilized from drugs
    clear_time: int
         number of months the drugs is effective or present in the host
    clear_count: int
         count of months since drug administered
    dfAdult: dataframe
         pandas dataframe containing adult parasites
    dfMF: dataframe
          pandas dataframe containing larval stage parasites     
         
    Returns
    ------
    dfAdult
    dfMF
    '''
    if clear_count == 1: #permanent sterility
        sterile = np.random.random(1,len(dfAdult.age))
        #here dfAdult.selF is the sum of allele effect on phenotype
        # 1 is wildtype, <1 is less fit than wt, >1 is more fit than wt
        if sterile <= sterile_p ** (dfAdult.selF):
             #permanently sterilized adults
             dfAdult.fec[[sterile]] = 0
             
    elif clear_count < clear_time: #Drugs cause temporary sterility over clear_time
        #linear function defining fecundity during drug clearance
        m = float(fecund - 1) / (clear_time - 1 )
        b = 1 - m * 1
        #new base fecundity under drugs
        sterile_t = (m * clear_count + b)
        #assign value to dfAdult.fec
        ageltsixIDX = which(dfAdult.age < 6 && dfAdult.fec != 0)
        dfAdult[ageltsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson((sterile_t ** dfAdult.selF)))

        #assign value to dfAdult.fec with older worms
        agegtsixIDX = which(dfAdult.age >= 6  && dfAdult.fec != 0)
        #linear function of reproductive decline         
        mage = float(0 - fecund) / (21 - 6)
        bage = 0 - m * 21
        #sterile_t = ((float((mage * dfAdult.age[agegtsix] + bage)-1)/(clear_time-1)) 
         #    * clear_count + (1 - (float((mage * dfAdult.age[agegtsix] + bage)-1) 
          #        / (clear_time-1))))
        dfAdult[agegtsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(
             (((float((mage * dfAdult.age[agegtsix] + bage)-1)/(clear_time-1)) 
             * clear_count + (1 - (float((mage * dfAdult.age[agegtsix] + bage)-1) 
                  / (clear_time-1)))) ** dfAdult.selF)))
        
    else: #base fecundity when no drugs   
        #find all indexes where age is less than 6 and adult is not premanently sterilized
        ageltsixIDX = which(dfAdult.age < 6 && dfAdult.fec != 0)
        avg_fitness = mean(dfAdult.selF)  #average population fitness
        #assign value to dfAdult.fec
        dfAdult[ageltsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson(fecund ** (1 - abs(avg_fitness - dfAdult.selF))))
        #find all indexes (rows) where age is greater than equal to 6 and not sterilized
        agegtsixIDX = which(dfAdult.age >= 6  && dfAdult.fec != 0)
        #linear function of reproductive decline         
        m = float(0 - fecund) / (21 - 6)
        b = 0 - m * 21
        #assign value to dfAdult.fec from linear function     
        dfAdult[agegtsixIDX] = dfAdult.assign(fec=lambda dfAdult: np.random.poisson((m 
              * dfAdult.age[agegtsixIDX] + b) ** (1 - abs(avg_fitness - dfAdult.selF))))
    
    #actual reproduction just copies adult info from dfAdult to dfMF
    #number of times indicated in the dfAdult.fec column
    dfMF = dfMF.append(rep(dfAdult,dfAdult.fec))

    return dfAdult, dfMF
 

if __name__ == '__main__':
    #fecundity_basefx(fecund=20, dfAdult, dfMF)
    #fecundity_mdafx(fecund=20, sterile_p=0.35, clear_time=6, clear_count=1, dfAdult, dfMF)
    #fecundity_mdasel1fx(fecund=20, sterile_p=0.35, clear_time=6, clear_count=1, dfAdult, dfMF)
    #fecundity_mdasel2fx(fecund=20, sterile_p=0.35, clear_time=6, clear_count=1, dfAdult, dfMF)
    

