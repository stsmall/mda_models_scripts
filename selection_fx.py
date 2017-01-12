# -*- coding: utf-8 -*-

import numpy as np

def dfe_fx(gamma_alpha=1, gamma_beta=4):
     '''function to determine the fitness effect of a new mutation
     Parameters
     --------
     gamma_alpha
     gamma_beta
     
     Returns
     ------
     fitness value as float
     '''
     
     return np.random.gamma(gamma_alpha,gamma_beta)
     