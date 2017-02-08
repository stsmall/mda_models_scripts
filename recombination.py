#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
FiGS Copyright (C) 2017 Scott T. Small
This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
This is free software, and you are welcome to redistribute it
under certain conditions; type `show c' for details.
'''
import numpy as np
import random
import pandas as pd
from IPython import embed

def recombination_locus(h1, h2, basepairs):
    """Calculates the recombination at a given locus
    """
    crossover_pos = random.randint(0, basepairs)
    h1_ix = [i for i, x in enumerate(h1) if x > crossover_pos][0]
    h2_ix = [i for i, x in enumerate(h2) if x > crossover_pos][0]
    h1_new = np.append(h1[0:h1_ix + 1], h2[h2_ix:])
    h2_new = np.append(h2[0:h2_ix + 1], h1[h1_ix:])
    return(h1_new, h2_new)


def recombination_fx(locus,
                     dfAdult,
                     recombination_rate,
                     basepairs):
    """calculate number of recombination events and rearranges haplotypes

    Parameters
    ---------
    locus: int
        number of loci
    dfAdult_mf : pandas dataframe
        dataframe containing new larval parasites
    dfAdult : pd df
        dataframe containing reproducing adults
    recombination_rate : float, list
        recombination rate for each locus
    basepairs : int, list
        length of each locus in basepairs

    Returns
    -------
    dfAdult_mf : pd df

    """
    males = dfAdult.query('sex == "M"')
    newmf = []
    # For each female
    for index, row in dfAdult.query('sex == "F"').iterrows():
        hostidx = row.hostidx
        try:
            male = males.query('hostidx == @hostidx').sample(1)
        except ValueError:
            print("No males in human host, no sex")
            continue
        new_row = row.copy()
        for mf in range(row.fec):
            for loc in range(locus):
                # locus identifier
                lid = "locus_" + str(loc) + '_h{0!s}'
                if recombination_rate[loc] == 0:
                    continue
                else:
                    num_recomb = np.random.poisson(
                            recombination_rate[loc] * basepairs[loc] * 2)
                    if num_recomb == 0:
                        new_row[lid.format(1)] =\
                                row[lid.format(random.choice("12"))]
                        #male contribution
                        new_row[lid.format(2)] =\
                                male[lid.format(random.choice("12"))].values[0]
                    else:
                        #randomly choose male or female
                        sex_xing = random.choice("MF")
                        # :TODO need to fix how this is initialized
                        h1m = male[lid.format(1)].values[0].copy()
                        h2m = male[lid.format(2)].values[0].copy()
                        h1f = row[lid.format(1)].copy()
                        h2f = row[lid.format(2)].copy()
                        if sex_xing is "M":
                            for _ in range(num_recomb):
                                h1m, h2m = recombination_locus(h1m, h2m,
                                        basepairs[loc])
                        elif sex_xing is "F":
                            for _ in range(num_recomb):
                                h1f, h2f = recombination_locus(h1f, h2f,
                                        basepairs[loc])
                        new_row[lid.format(1)] = random.choice([h1f, h2f])
                        new_row[lid.format(2)] = random.choice([h1m, h2m])
                    newmf.append(new_row)
    dfAdult_mf = pd.concat(newmf, ignore_index=True, axis=1).T
    return dfAdult_mf
