# Wb simulation model
Simulation to generate genetic sequence data for specific life history of Wuchereria bancrofti.

# Documentation
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

### Function tree
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
