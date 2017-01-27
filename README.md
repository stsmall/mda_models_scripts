# Filarial worm genetic simulation model *(FiGS)*
There are currently no available simulation models that generate population genetic data for filarial worms. Population genetic data can be used to infer specific life history characterisitcs of filarial worm populations and track the response and progress of elimination efforts. FiGS is a forward-in-time simulation written in python to generate genetic data under the specfic life history of *W. bancrofti*. FiGS is currently being used to identify sutaible population genetic summary statistics for evaluating elimination progression. Although written specifically for *W. bancrofti*, it is probably tunable to other filarial worms by altering the historical demographic parameters in the coalescent functions.

**Relevant publications**:   
[Ramesh A, *et al.,* **(2012)**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3725818/): 1st published mitochondrial genome for *W. bancrofti*  
[Small, ST, *et al.,* **(2013)**](http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0002308): population genetics of *W. bancrofti* from Papua New Guinea  
[Small, ST, *et al.,* **(2014)**](http://europepmc.org/articles/pmc4257870): review of techniques for studying genetic diversity in *W. bancrofti*  
[Small, ST, *et al.,* **(2016)**](http://onlinelibrary.wiley.com/doi/10.1111/mec.13574/full): 1st published population genomics in *W. bancrofti*  

**Current Assembly**  
[wormbase](http://parasite.wormbase.org/Wuchereria_bancrofti_prjna275548/Info/Index): Whole Genome Assembly from Small 2016     

# Documentation
## Dependencies
Anaconda will install all dependencies as well as scrm (e.g., conda install -c bioconda scrm)  
* [scrm](https://scrm.github.io/)
* scipy, numpy, sklearn, pandas, matplotlib
* [scikit-allel](https://scikit-allel.readthedocs.io/en/latest/)

## options
full options list available in docstrings

#### Host Demography
* villages: number of villages to simulate
* host population size: population of each village
* prevalence: initial prevalence of W.bancrofti infection
* host migration rates: movement of hosts between villages

#### Vector
* biting rates: bites per person per hour
* biting times: avg exposure to biting

#### Genetic
* number of loci: how many loci to simulate
* length of loci: how many basepairs. scrm can handle Mb size loci with recombination
* mutation and recombination rates: locus specific
* selection:

#### Intervention
* mass drug administration: MDA
* bed nets: 

#### Data recording
* demography: formatted to R dataframe
* population genetic: summary statistics from pylibseq
* genetic data: output all or sample as vcf

# Functions tree
#### wbsims_initialize.py
* wbsims_init
 * agehost_fx
 * host_fx
 * coalsims_migmatrix_fx
 * parse_coalsims_fx
 * coalsims_fx
 * sel_fx
 * fit_fx
 * wormdf_fx

#### wbsims_run.py
* wb_sims
  * transmission.py
    * hostmigration_fx
    * vectorbite_fx
    * agehost_fx
    * newinfection_fx
    * transmission_fx
  * survival.py
    * survivalbase_fx
    * fecundity.py
      * fecunditybase_fx
        * recombination.py
          * recombination_fx
        * mutation.py
          * mutation_fx
          * selection.py
            * selection_fx
            * fitness_fx
  * record_data.py
    * record_demo_fx
    * record_popgen_fx
      * gt2scikit_fx
        * scikit-allel
# License        
    Filarial worm genetic simulations, FiGS
    Copyright (C) 2017  Scott T. Small

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see [license](http://www.gnu.org/licenses/)
