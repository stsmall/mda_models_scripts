# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 18:51:45 2017

@author: stsmall
"""
###conditions mda=T, selection=T    

import numpy as np
import pandas
from scipy.stats import weibull_min

dfAdult = pd.DataFrame(
            {"village" : [],
             "host" : [],
             "age" : [],
             "R0net" : [],
             "fec" : [],      #records fecundity, lambda for poisson()
             "hap" : [],
             "dip1" : [],
             "dip2" : [],
             "selF" : [],    #records fitness offset on fecundity
             "selS" : []     #records fitness offset on survival
             })
dfJuv = pd.DataFrame(
            {"village" : [],
             "host" : [],
             "age" : [],
             "R0net" : [],             
             "hap" : [],
             "dip1" : [],
             "dip2" : [],
             "selF" : [],    #records fitness offset on fecundity
             "selS" : []     #records fitness offset on survival
             })
dfMF = pd.DataFrame(
            {"village" : [],
             "host" : [],
             "age" : [],
             "R0net" : [],             
             "hap" : [],
             "dip1" : [],
             "dip2" : [],
             "selF" : [],    #records fitness offset on fecundity
             "selS" : []     #records fitness offset on survival
             })
             
dfSel = pd.DataFrame(
            {"locus" : [],
             "position" : [],
             "selF" : [],
             "selS" : []})

def fecundity_mdaselfx(fecund=20, 
                    sterile_p=0.35, 
                    clear_time=6, 
                    clear_count=1):

if __name__ == '__main__':
    fecundity_mdafx(fecund=20, sterile_p=0.35, clear_time=6, clear_count=1)
