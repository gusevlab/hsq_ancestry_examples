# Heritability analyses

`reml_func.R` : Utility functions for running GCTA-style REML and HE regression in R

`sib_RDR_comparison.R` : Script to simulate family data and estimate direct h2 using RDR and a sibling version of RDR.

Running without any parameter changes will simulate five instances of populations with direct h2 of 0.2, indirect variance of 0.1, and then estimate the direct h2 with different methods. The output looks like this:

```
Mean_h2	Mean_SE	SE_h2
True	0.20321946241417	NA	0.00494282242025141
POP_h2	0.451991484291674	NA	0.00799286062770628
RDR_fam	0.186236463907084	0.0342745893623424	0.0148506914090754
RDR_sib	0.183529708989363	0.0543929865454116	0.0216861408233236
```

Where the first column is the method type, the second column is the mean direct h2 estimate (or the population estimate for the "POP_h2" row), the third column is the average standard error computed by GREML, and the fourth column is the standard error of the mean direct h2 estimate over the simulated instances.