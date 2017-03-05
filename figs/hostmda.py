#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np
import ipdb

def hostmda_fx(village,
               dfHost,
               mda_coverage):
     ''' Assigns MDA to hosts per village based on mda_coverage
     Parameters
     ----------
     villages : int
          number of villages
     dfHost : df
          dataframe containing host information
     mda_coverage : list, float
          coverage of mda for each village

     Returns
     -------
     dfHost : df
          now will mda column as identity 0 or 1, and updates mda_sum column
     '''
     #ipdb.set_trace()
     y = lambda z: z.append(1)
     x = lambda z: z.append(0)
     ##add to MDA_cum
     dfHost.loc[dfHost.MDAstate == 1].MDAcum.apply(y)
     dfHost.loc[dfHost.MDAstate == 0].MDAcum.apply(x)
     ##rezero MDA identifier
     dfHost.MDA = 0
     ##assign MDA randomly
     for vill in range(len(village)):
         mdaiix = dfHost[dfHost.village == vill].sample(frac = mda_coverage[vill]).index.values
         dfHost.loc[mdaiix, "MDAstate"] = 1
     return(dfHost)
