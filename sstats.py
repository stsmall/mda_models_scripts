#!/usr/bin/env python
from pandas import *
import sys
def stats(file):
    stats_ms=read_table(file)
    mean=stats_ms.mean()
    var=stats_ms.var()
    return mean
print stats(sys.stdin)