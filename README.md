# Target capture enrichment sequencing of SNP loci

This repository contains the custom scripts for the paper


Pacheco C, Diana L, Silva P, Álvares F, García EJ, Castro D, Layna J, López-Bao JV & Godinho R (2022) Is there a tissue of choice of museomics? *In prep*


The scripts are:

* gtvalues2plink.py to create a PLINK dataset (.map & .ped files) from the tabular output of GATK VariantToTable, with the option of filtering based on minimum values;

* SNP_concordance.py to calculate the concordance rate between two PLINK datasets.


See the comments within each script for instructions about their usage.