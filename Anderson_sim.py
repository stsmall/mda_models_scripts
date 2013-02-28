import numpy
import math
import stat

transmission=open("transmission_dist.txt", "a")
dna=open("outfile","r")
distances=[]
within=[]
among=[]
for line in dna:
    for item in line.split():
        if "." in item:   #index
            distances.append(float(item))
for i in range(1,40):
    temp=distances[(40*i):(40*(i+1))]
    if i <20:
        within.extend(temp[0:19])
        among.extend(temp[20:40])
    else:
        within.extend(temp[20:40])
        among.extend(temp[0:19])
    
#append this to a file
transmission.write("\n%s" %(numpy.mean(within)/numpy.mean(among)))
transmission.close()
#within=Series(within)
#among=Series(among)
#among.to_csv('Anderson_among.txt')
#within.to_csv('Anderson_within.txt')

#np.asarray(within).T  #transpose
