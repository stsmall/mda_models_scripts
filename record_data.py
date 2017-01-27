#!/usr/bin/env python2
# -*- coding: utf-8 -*-
'''
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
'''

def allelefreq_fx():
          allele = []
    
     ##freq of positions in dfSel
     for loc in range(1, locus):
          for row in range(len(dfAdult)):
               #flatten list and extend
               allele.extend(dfAdult.iloc[row]["locus_" + str(loc) + "_h1"])
               allele.extend(dfAdult.iloc[row]["locus_" + str(loc) + "_h2"])
     freq = np.unique(allele, return_counts=True)[1] / (len(dfAdult) * 2.0) 