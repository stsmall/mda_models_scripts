# Wb simulation model
Simulation to generate genetic sequence data for specific life history of Wuchereria bancrofti.

# Documentation
## options [^1]

#### Host Demography
##### villages
##### host population size
##### prevalence
##### host migration rates

#### Vector
##### biting rates
##### biting times

#### Genetic
##### number of loci
##### length of loci
##### mutation and recombination rates
##### selection

#### Intervention
##### mass drug administration MDA
##### bed nets

[^1]: full options list availabe in docstrings

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
