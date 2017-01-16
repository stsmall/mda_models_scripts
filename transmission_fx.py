# -*- coding: utf-8 -*-
import numpy as np
from sklearn.metrics import pairwise_distances
from random import randrange  
from math import sqrt
def transmission_fx(L3trans=100, sigma=50, dfMF, dfJuv, dfHost):
    '''transmission events resolved as either reinfection or new infection
    Parameters
    --------
    L3trans:int 
         number of successful transmission events from vector fx
    dispersal: float
         dispersal distance of mosquito, circle with radius 2*sigma
    dfHost: df
         chooses the donating and accepting hosts for transmission
    dfMF: df
         donating MF genotype to transmit
    dfJuv: df
         donated MF move to Juvenille age class

    Returns
    ------
    dfJuv
    
    '''
    
    #choose host donating MF, weighed by number of MF
    donatingMF = np.random.randint(0,len(dfMF.age), L3trans)  
    #distance between pair of hosts as sparse matrix
    #dist = math.sqrt((x2[0] - x1[0]) ** 2 + (y2[1] - y1[1]) ** 2)
    distMat = pairwise_distance([dfHost.xpos,dfHost.ypos])
    #determines how far mosquito can transmit
    dispersal = 2 * sigma    
    
    #if new_infection w/ prob
    if len(dfHost.age) < hostpopsize:
         prob_newinfection = float(1)/np.sum(distMat[dfHost.where(dfMF.hostidx[dontaingMF])] < dispersal)
    else:
         prob_newinfection = 0  #maintains zero-sum
    
    if np.random.random() < prob_newinfection:     
         #create new host, drawing values for age, agedeath, location. where location
         dfHost, newhostidx = new_infection_fx(mindist=1, dispersal, dfMF.idx[dontatingMF], dfHost)
         mftrans = dfMF[donatingMF]
         mftrans.hostidx = newhostidx
         dfJuv.append(mvtrans)
    #reinfection
    else:
         #accepting host <= dispersal distance from donating host
         acceptingHost = np.random.randint(0,len(dfHost.age))
         #add row to dfJuv, with idx matching new host
         mftrans = dfMF[donatingMF]
         mftrans.hostidx = dfHost.hostidx[acceptingHost]
         dfJuv.append(mftrans)    
         
    return dfHost, dfJuv

def new_infection_fx(mindist=1, dispersal, donatingidx, dfHost):
    '''a transmission event infecting a naive host

    Parameters
    ----------
    dispersal: float
         dispersal distance as 2*sigma
    dfHost: df
        host dataframe 
    donatingidx: str
        index of host where MF is originating

    Returns
    ------
    dfHost: df
         updated with new host information
    newhostidx: str
          new host index     
    '''
    ##assign new village
     #copy from donating host
    ##assign new hostidx
    
    
    ##assign new sex
    random.choice("MF")
    ##assign new age
    #see intial functions
    
    ##assign new death
    #see intial functions
    
    #assign new position
    maxc=int(sqrt((maxdist**2/2)))
    x_new = randrange(-maxc,maxc,1)
    y_new = randrange(-maxc,maxc,1)
    if sqrt(x_new**2+y_new**2) > mindist:
         hostptX = dfHost.pos[dfHost.hostidx[donatingidx]][0] + x_new
         hostptY = dfhost.pos[dfHost.hostidx[dontaingidx]][1] + y_new
         newpts = [hostptX, hostptY]
         if newpts not in dfHost.pos[:]:
              newhost = dfHost.hostidx[donatingidx]
              newhost[1]
              dfHost.append(newhost)
    
    return dfHost, newhostidx
    