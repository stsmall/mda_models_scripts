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
     if move_host != 0:
          migrant = np.random.randint(0, len(dfHost), move_host)
          for mv in migrant:
               if dfHost.loc[migrant, "village"] < max(dfHost.village):
                    dfHost.loc[migrant, "village"] += 1
               else:
                    dfHost.loc[migrant, "village"] -= 1
     return(dfHost)
