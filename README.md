# Wb simulation models
Simulation to generate genetic sequence data for specific life history of Wuchereria bancrofti.

##Documentation
###options

###functions
####wbsims_initialize.py
* wbsims_init_fx
 * host_fx
  * agehost_fx
 * wormdf_fx
  * ms_outcall
    * migration_matrix
    * parse_ms_output
  * sel_fx
  * fitness_fx

####wbsims_run.py
* maturation_fx
  * survival_fx
  * fecundity_fx
    * mutation_fx
    * recombination_fx
  * sel_fx
  * fitness_fx   
* vectorbite_fx
  * transmission_fx
* write_demog_fx
* write_popgen_fx
  * gt2scikit_fx
