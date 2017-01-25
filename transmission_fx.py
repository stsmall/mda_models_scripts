# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import math
from sklearn.metrics import pairwise_distances
import random

def vectorbite_fx(bitespperson, 
                  hours2bite, 
                  hostpopsize,
                  prev_t,
                  densitydep, 
                  avgMF, 
                  bednets, 
                  bnstart, 
                  bnstop, 
                  bncoverage, 
                  month):
    '''counts the number of successful infectious mosquito bites
    Parameters
    --------
    bitespperson : int
        rate of bites per person per unit time, here hours
    hours2bite : int
        number of hours mosq bite per night/day
    hostpopsize : int
        total population of village 
    prev_t : float
        percent of people infected 
    densityDep : Boolean
        use the density dependece function for developing L3
    avgMF : float
         average number of MF per host per village
    bednets : Boolean, list
        use values from bednets     
    bnstart : int, list
        month to start bednets 
    bnstop : int, list
        mont to stop bednets  
    bncoverage : int, list
        percent of population using bednets 
    month : int
        current time, month 
    Returns
    ------
    L3trans : int
         all new juv are age class 0
    '''
    
    if bednets:
        if month > bnstart and month < bnstop:
             totalbites = ((1 - bncoverage) * bitespperson * hours2bite * 30 
                                   * hostpopsize)
        else:
             totalbites = (bitespperson * hours2bite * 30 
                                   * hostpopsize)
    else:
        totalbites = (bitespperson * hours2bite * 30 
                              * hostpopsize)
    # 0.37 is prob of bite on infected host picking up MF    
    infbites = np.random.binomial(totalbites, (prev_t * 0.37))
    if densitydep: #values for anopheles from CITE
       #number of MF in 20ul of blood
       #235ml is 5% of total host blood 
       #there are 50 units of 20ul in 1ml
       mfBlood = float(avgMF/ 235) / 50 
       # 0.414 is proportion of  L3 that leave mosquito per bite 
       # 0.32 proportion of L3 that make it into the host
       L3trans = round(infbites * (4.395 * (1 - math.exp( -(0.055 * (mfBlood)) 
           / 4.395)) ** 2) * (0.414 * 0.32)) 
    else:
       # 0.414 is proportion of  L3 that leave mosquito per bite 
       # 0.32 proportion of L3 that make it into the host
       L3trans = np.random.binomial(infbites, (0.414 * 0.32))
    
    return L3trans

def new_infection_fx(mindist, dispersal, dfHost, transMF):
    '''a transmission event infecting a naive host

    Parameters
    ----------
    mindist : int
         minimum distance
    dispersal: float
         dispersal distance as 2*sigma
    dfHost: df
        host dataframe 
    transMF
        row and all columns of transmitted MF
    Returns
    ------
    dfHost: df
         updated with new host information
    newhostidx: str
          new host index     
    '''
    #assign new position
    maxc = int( math.sqrt(( maxdist ** 2 / 2 )))
    x_new = random.randrange( -maxc, maxc, 1 )
    y_new = random.randrange( -maxc, maxc, 1 )
    if sqrt(x_new ** 2 + y_new ** 2) > mindist:
         hostptX = dfHost[dfHost.hostidx == transMF.hostidx].coordinates[0] + x_new
         hostptY = dfHost[dfHost.hostidx == transMF.hostidx].coordinates[1] + y_new
         newpts = [hostptX, hostptY]
         #check that new position does not already exist
         if newpts not in dfHost.coordinates[:]:
              #copy village
              vill = transMF.village
              #new host index
              oldhostidx = dfHost[dfHost.village == transMF.village].hostidx.iloc[-1]
              new_hostidx = oldhostidx[:oldhostidx.rfind('h')] + 'h' + str(int(oldhostidx.split('h')[1]) + 1)
              #radnom sex
              sex = random.choice("01")
              #age and agedeath function from wbsims_initialize
              age, agedeath = agehost_fx(sex)
              #assign coordinates
              newpts
              #add to df
              dfHost.loc[-1] = [vill, new_hostidx, sex, age, agedeath, newpts, 0]
              dfHost.index = dfHost.index + 1
              dfHost = dfHost.sort()
              
    return dfHost, new_hostidx

def transmission_fx(villages, hostpopsize, sigma, bitesPperson, hours2bite, densityDep, bednets,
                    bnstart, bnstop, bncoverage, month, dfMF, dfJuv, dfHost):
     #transmission_fx(2, [200, 200], 50, [20, 20], [8, 8], [True, True], [False, False], 
                      #[12, 12], [36, 36], 1, dfMF, dfJuv, dfHost) 
    '''transmission events resolved as either reinfection or new infection
    Parameters
    --------
    villages : int
        number of villages
    hostpopsize : int, list
        total population size of village      
    sigma : float
        dispersal mean
    bitespperson : int, list
        rate of bites per person per unit time, here hours
    hours2bite : int, list
        number of hours mosq bite per night/day
    densityDep : Boolean
        use the density dependece function for developing L3
    bednets : Boolean, list
        use values from bednets     
    bnstart : int, list
        month to start bednets 
    bnstop : int, list
        mont to stop bednets  
    bncoverage : int, list
        percent of population using bednets 
    month : int
        current time, month 
    dfHost: df
         chooses the donating and accepting hosts for transmission
    dfMF: df
         donating MF genotype to transmit
    dfJuv: df
         donated MF move to Juvenille age class
    Returns
    ------
    dfJuv : df
         all new juv are age class 0
    dfHost : df
         in the case of newinfection, new host
    '''
    dispersal = 2 * sigma
    for vill in range(villages):
         prev_t = len(dfHost[dfHost.village == vill])/hostpopsize[vill]
         infhost = len(dfHost[dfHost.village == vill])
         avgMF = len(dfMF[dfMF.village == vill])/float(infhost)
         L3trans = vectorbite_fx(bitesPperson[vill], hours2bite[vill], hostpopsize[vill], 
                                 prev_t, densityDep[vill], avgMF, bednets[vill], 
                                 bnstart[vill], bnstop[vill], bncoverage[vill], month)    
         if L3trans > len(dfMF[dfMF.village == vill]):
              transMF = dfMF[dfMF.village == vill]     
         else:
              transMF = dfMF[dfMF.village == vill].sample(L3trans)
              
         distMat = pairwise_distances(np.vstack(dfHost[dfHost.village == vill].coordinates))                  
         
         for index, row in transMF.iterrows():
              #new infection
              if len(dfHost[dfHost.village == vill]) < hostpopsize[vill]:
                   prob_newinfection = 1.0/(len(distMat[index] < dispersal) - 1)
              else:
                   prob_newinfection = 0
              
              if np.random.random() < prob_newinfection:     
                   #new host
                   dfHost, newidx = new_infection_fx(dispersal, row, dfHost)
                   row.hostidx = newidx
                   row.age = 0
                   dfJuv.append(row)
              else: #reinfection
                   rehost = dfHost.loc[random.choice(np.where(distMat[index] <= dispersal)[0])]
                   row.hostidx = rehost.hostidx
                   row.age = 0
                   dfJuv.append(row)
                            
    return dfHost, dfJuv