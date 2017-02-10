#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 16:33:14 2017

@author: scott
"""

def counthaps(hap):
    a=[]
    [a.append(i) for i in hap]
    b=[]
    [b.append(str(list(i))) for i in a]
    counts = dict()
    for i in b:
        counts[i] = counts.get(i,0) + 1
    return(counts.values())