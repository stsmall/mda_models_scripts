# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 12:19:29 2017
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

###conditions: mda=False, selection=False
    
def survival_basefx(month=1, 
                    surv_Juv=0.866, 
                    shapeMF=3.3, 
                    scaleMF=12, 
                    shapeAdult=3.8, 
                    scaleAdult=8):
                        
    '''base survival function'''
    #conditions: mda=False, selection=False
    ##juv is exponential 0.866   
    survjuv = np.random.random(len(dfJuv.age)) #array of random numbers
    killjuv = which(survjuv >= surv_Juv)    
    dfJuv.drop(dfJuv[[killjuv]])
    
    ##MF is weibull cdf
    survmf = np.random.random(len(dfMF.age))    
    surv_MF = apply(dfMF.age, 2, function(x) weibull_min.cdf(dfMF.age,shapeMF,
                    loc=0,scale=scaleMF))      
    killmf = which(survmf <= surv_MF)     
    dfMF.drop(dfMF[[killmf]]) 

    if month%12 == 0:   
        ##Adult is weibull cdf
        survadult = np.random.random(len(dfAdult.age))    
        surv_Adult = apply(dfAdult.age, 2, function(x) weibull_min.cdf(dfAdult.age,
                           shapeAdult,loc=0,scale=scaleAdult))       
        killadult = which(survadult <= surv_Adult)     
        dfAdult.drop(dfAdult[[killadult]]) 
        
if __name__ == '__main__':
    survival_basefx(surv_Juv=0.866, shapeMF=3.3, scaleMF=12, shapeAdult=3.8, scaleAdult=8)




    