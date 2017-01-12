# Wb simulation models
Simulation to generate genetic sequence data for specific life history of Wuchereria bancrofti.

##Documentation
###Parameters
####Host
* Villages
* MDA
* Bednet
* Prevalence
* hostpopsize

####Parasite


####Genetics
* selection



###Data Frame structures
####data frame for host rates
dfHost = pd.DataFrame(
            {"village" : [], #village identifier as int
             "host" : [], #host identifier as int within Villages
             "age" : [], #age of host for life table of morality
             "angle" : [], #polar coordinate for location
             "radius" : [], #polar cooredinate for location
             "MDA" : [] #participation in MDA, Yes = 1, No = 0
             })
####data frame for adult parasites
dfAdult = pd.DataFrame(
            {"village" : [], #village identifier as int
             "host" : [], #host identifier as int within villages
             "age" : [], #age of parasite in years
             "R0net" : [], #random number to keep track of net reproductive fitness
             "fec" : [], #rate parameter for poisson distribution
             "loc1" : [], #genotype of haploid locus, e.g., [9,234]
             "loc2" : [], #genotype of diploid locus, e.g., [[9.234],[10,456]]
             "loc3" : [], #genotype of diploid locus, e.g., [[9.234],[10,456]]
             "selF" : [], #net fitness for alleles contributing to fecundity
             "selS" : [] #net fitness for alleles contributing to survival
             })
####data frame for juvenille parasites
dfJuv = pd.DataFrame(
            {"village" : [], #village identifier as int
             "host" : [], #host identifier as int within villages
             "age" : [], #age of parasite in months
             "R0net" : [], #random number to keep track of net reproductive fitness         
             "hap" : [], #genotype of haploid locus, e.g., [9,234]
             "dip1" : [], #genotype of diploid locus, e.g., [[9.234],[10,456]]
             "dip2" : [], #genotype of diploid locus, e.g., [[9.234],[10,456]]
             "selF" : [], #net fitness for alleles contributing to fecundity
             "selS" : [] #net fitness for alleles contributing to survival
             })
####data frame for larval parasites; microfilaria
dfMF = pd.DataFrame(
            {"village" : [], #village identifier as int
             "host" : [], #host identifier as int within villages
             "age" : [], #age of parasite in months
             "R0net" : [], #random number to keep track of net reproductive fitness             
             "hap" : [], #genotype of haploid locus, e.g., [9,234]
             "dip1" : [], #genotype of diploid locus, e.g., [[9.234],[10,456]]
             "dip2" : [], #genotype of diploid locus, e.g., [[9.234],[10,456]]
             "selF" : [], #net fitness for alleles contributing to fecundity
             "selS" : [] #net fitness for alleles contributing to survival
             })
####data frame for tracking phenotype of mutations
dfSel = pd.DataFrame(
            {"locus" : [], #locus number
             "position" : [], #position of mutation, allele
             "selF" : [], #fitness effect on fecundity
             "selS" : []}) #fitness effect on survival

###Function descriptions
####
