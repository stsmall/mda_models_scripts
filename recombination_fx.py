import random


def recombination_fx(mf, num_recomb, basepairs):
    """this is run every time the prob of a recombination is true

    Parameters
    ----------
    mf: list
        list of mf from adults
    num_recomb: list of ints
        number of recombination events observed
    basepairs:

    Returns
    -------
    """
    for i, bp in enumerate(basepairs):
        recomb = 0
        while recomb < num_recomb[i]:  
            # choose random index in mf lis
            rec_mf = random.randrange(0, len(mf))
            new_recomb = random.randint(0, 3)
            if new_recomb < 2:  # first parent
                hap1 = mf[rec_mf][i + 1][0]
                hap2 = mf[rec_mf][i + 1][1]
                crossover_pos = random.randint(0, bp)
                hap1.sort()
                hap2.sort()
                hap1_co = next(l[0] for l in enumerate(
                    hap1) if l[1] > crossover_pos)
                hap2_co = next(l[0] for l in enumerate(
                    hap2) if l[1] > crossover_pos)
                hap1_new = hap1[0:hap1_co] + hap2[hap2_co:]
                hap2_new = hap2[0:hap2_co] + hap1[hap1_co:]
                mf[rec_mf][i + 1][0] = hap1_new
                mf[rec_mf][i + 1][1] = hap2_new
            else:  # second parent
                hap3 = mf[rec_mf][i + 1][2]
                hap4 = mf[rec_mf][i + 1][3]
                crossover_pos = random.randint(0, bp)
                hap3.sort()
                hap4.sort()
                hap3_co = next(l[0] for l in enumerate(
                    hap3) if l[1] > crossover_pos)
                hap4_co = next(l[0] for l in enumerate(
                    hap4) if l[1] > crossover_pos)
                hap3_new = hap1[0:hap3_co] + hap2[hap4_co:]
                hap4_new = hap2[0:hap4_co] + hap1[hap3_co:]
                mf[rec_mf][i + 1][2] = hap3_new
                mf[rec_mf][i + 1][3] = hap4_new
            recomb += 1
    return mf
