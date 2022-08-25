# Scripts for Pacheco, Lobo et al. 2022: "Assessing the performance of historical skins and bones for museomics using wolf specimens as a case study"

This repository contains the custom scripts for the paper

```
Pacheco C, Lobo D, Silva P, Álvares F, García EJ, Castro D, Layna JF, López-Bao JV and Godinho R (2022)
Assessing the performance of historical skins and bones for museomics using wolf specimens as a case study.
Frontiers in Ecology Evolution 10:970249
```
[DOI: 10.3389/fevo.2022.970249](https://doi.org/10.3389/fevo.2022.970249)

The scripts are:

* `gtvalues2plink.py` to create a PLINK dataset (.map & .ped files) from the tabular output of GATK VariantToTable, with the option of filtering based on minimum values;

* `SNP_concordance.py` to calculate the concordance rate between two PLINK datasets.


See the comments within each script for instructions about their usage.
