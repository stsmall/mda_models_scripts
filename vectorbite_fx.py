# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 08:49:38 2017

@author: stsmall
"""
import numpy as np
import math
def vectorbite_fx(villages=2, 
                  bitespperson=[20, 20], 
                  hours2bite=[8, 8], 
                  hostpopsize=[100, 150], 
                  prev_t=[0.3, 0.2], 
                  densitydep=True, 
                  villavgMF=[65, 90]):
    '''counts the number of successful infectious mosquito bites

    Parameters
    ----------
    villages: int, list
        number of villages as grouping
    bitespperson: int, list
        rate of bites per person per unit time, here hours
    hostpopsize: int, list
        size of the host population in each village
    prev_t: float, list
        the prevalence at time t (per month)
    densitydep: Boolean
        use the density dependece function for developing L3
    villavgMF: float, list
        avgerage number of MF per host per village
    
    Returns
    ------
    transL3
    
    ---'''
                  
    for rate in range(villages):
        # 30 days per month    
        totalbites_village = (bitespperson[rate] * hours2bite[rate] * 30
            * hostpopsize[rate])
        # 0.37 is prob of bite on infected host picking up MF    
        infbites = np.random.binomial(totalbites_village, (prev_t[rate] * 0.37))
        if densitydep: #values for anopheles from CITE
            #number of MF in 20ul of blood
            #235ml is 5% of total host blood 
            #there are 50 units of 20ul in 1ml
            mfBlood = float(villavgMF[rate] / 235) / 50 
            # 0.414 is proportion of  L3 that leave mosquito per bite 
            # 0.32 proportion of L3 that make it into the host
            transL3 = round(infbites * (4.395 * (1 - math.exp( -(0.055 * (mfBlood)) 
                / 4.395)) ** 2) * (0.414 * 0.32)) 
        else:
            # 0.414 is proportion of  L3 that leave mosquito per bite 
            # 0.32 proportion of L3 that make it into the host
            transL3 = np.random.binomial(infbites, (0.414 * 0.32))
    
    return transL3
