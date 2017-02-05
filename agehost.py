#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np

def agehost_fx(sex):
    '''Calculate age of host and age of death from acturial table
    Parameters
    ---------
    sex : str
         male or female

    Returns
    ------
    age : int
         age of host
    death : int
         age of death
    '''
    zero_10 = 0.24 #12%
    eleven_20 = 0.46 #11%
    twentyone_30 = 0.64 #9%
    thirtyone_40 = 0.78 #7%
    fourtyone_50 = 0.88 #5%
    fiftyone_60 = 0.94 #3%
    sixtyone_70 = 0.98 #2%                     
    #assign age
    agerand = np.random.random()
    if agerand <= zero_10:
         age = np.random.randint(1,10)
    elif agerand > zero_10 and agerand <= eleven_20:
         age = np.random.randint(11,20)
    elif agerand > eleven_20 and agerand <= twentyone_30:
         age = np.random.randint(21,30)
    elif agerand > twentyone_30 and agerand <= thirtyone_40:
         age = np.random.randint(31,40)
    elif agerand > thirtyone_40 and agerand <= fourtyone_50:
         age = np.random.randint(41,50)
    elif agerand > fourtyone_50 and agerand <= fiftyone_60:
         age = np.random.randint(51,60)
    elif agerand > fiftyone_60 and agerand <= sixtyone_70:
         age = np.random.randint(61,70)
    elif agerand > sixtyone_70:
         age = np.random.randint(71,80)
    
    #dictionary from actuarial tables 
    deathdict = {}
    with open("./act.tbl",'r') as tbl:
         for line in tbl:
              line = line.strip()
              deathdict["{}".format(line.split()[0])] = list(map(float,
                  line.split()[1:]))
    
    #when do they die
    death = deathdict[str(age)][int(sex)] + np.random.normal(0,6) + age
         
    return(age, round(death))
