#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np


def hostmigration_fx(dfHost,
                     hostmigrate):
     ''' allows host to move between villages
     Parameters
     ---------
     dfHost : df
          df of hosts
     hostmigrate : float
          host migration rate between villages per year
     Returns
     --------
     dfHost : df
     '''
     move_host = np.random.poisson(hostmigrate)
     i = 0
     while i < move_host:
         migrant = np.random.choice(dfHost.index, move_host)
         #stepping stone
         for mv in migrant:
             if dfHost.ix[migrant, "village"] == min(dfHost.village):
                  #less than max can go up or down
                  dfHost.ix[migrant, "village"] += 1
             elif dfHost.ix[migrant, "village"] == max(dfHost.village):
                  dfHost.ix[migrant, "village"] -= 1
             else: #at max can only go down
                 dfHost.loc[migrant, "village"] += np.random.choice([1,-1])
         i += 1
     return(dfHost)
