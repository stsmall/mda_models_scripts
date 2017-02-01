"""
    FiGS Copyright (C) 2017 Scott T. Small
    This program comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.
"""
import numpy as np

def filtercoords_fx(locus, positions, basepairs, perc_locus, cds_length, intgen_length):
    ''' Filters coordinates to keep only those in defined cds
    Parameters
    ---------
    
    Returns
    -------
    
    
    '''
    cds_positions = []
    cds_coordinates = []
    #size parameter
    size = 3
    #average distance between in bp
    mu = intgen_length 
    #last coordinate is approx num_cds * mu; so if num_cds is too large or mu is too long
    #genes will run over the locus length
    for loc in range(1,locus):
        num_cds = int(round((perc_locus[loc]*basepairs[loc]) / cds_length))
        size_cds = np.round(np.random.gamma(4, 0.25, num_cds) * cds_length)
        #r = size
        #m = mean
        #p = r / (  r + m )
        cds_between = np.random.negative_binomial(size, size/float(mu+size), num_cds)
        cds_stop = 0
        cds_coords = np.empty((0,2),int)
        for i, j in zip(cds_between, size_cds):
            #[i + cds_stop, i + cds_stop + j]
            cds_coords = np.append(cds_coords, [[i + cds_stop, i + j + cds_stop]], axis = 0)
            cds_stop += (i + j)
        #np.where(cds_coords > basepairs[loc])
        keep_pos = []
        for start, end in cds_coords:
            keep_pos.extend(positions[np.where(np.logical_and(positions >= start, 
                                                  positions <= end))]) 
        cds_positions.append(keep_pos)
        cds_coordinates.append(cds_coords)
    return(cds_positions, cds_coordinates)   