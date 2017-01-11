# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 19:23:35 2017

@author: stsmall
"""

def densitydep_fx(carrycap=100,):
                    Km = 100  #carrying capacity
                bm = 1.0/Km #Km is carrying capacity
                Sm = 0.99 #maximum survival
                am = (Sm*bm)/math.exp(-1)
                mort_A = am * sum_adult * math.exp(-bm * sum_adult) #Ricker fx
                mort_J = am * sum_juv * math.exp(-bm * sum_juv) #Ricker fx
                mort_M = am * sum_mf * math.exp(-bm * sum_mf) #Ricker fx
