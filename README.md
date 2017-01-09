
# PAINTOR
Probabilistic Annotation INtegraTOR

##UPDATE 01/07/17
Announcing PAINTOR v3.0. We have added new functionality to PAINTOR to conduct multi-trait fine-mapping. We have improved computational performance and numerical stability by employing a factorization of the MVN  In adddition, we have developed a visualiziation tool, PAINTOR-CANVIS, to produce publication-ready plots of an integrative fine-mapping experiment. 

## Description

We provide a command line implementation of the PAINTOR frameworks described in [Kichaev et al. (PLOS Genetics, 2014)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1004722)  [(American Journal of Human Genetics, 2015)](http://www.cell.com/ajhg/abstract/S0002-9297(15)00243-8).  [(Bioinformatics, 2016)](http://bioinformatics.oxfordjournals.org/content/early/2016/10/16/bioinformatics.btw615).  Briefly, PAINTOR is a statistical fine-mapping method that integrates functional genomic data with association strength from potentially multiple populations or traits to prioritize variants for follow-up analysis. The software runs on multiple fine-mapping loci and/or populations/traits simultaneously and takes as input the following data for each set of SNPs at a locus


1. Summary Association Statistics (Z-scores)
2. Linkage Disequilibrium Matrix/Matrices (Pairwise Pearson correlations coefficients between each SNP)
3. Functional Annotation Matrix (Binary indicator of annotation membership (i.e. if entry {i,k} = 1, then SNP i is a member of annotation K). 

#### Key Features

1. Outputs a probability for a SNP to be causal which can subsequently be used to prioritize variants
2. Can model multiple causal variants at any risk locus
3. Leverage functional genomic data as a prior probability to improve prioritization
  - This prior probability is not pre-specified, but rather, learned directly from the data via an Empirical Bayes approach
4. Quantify enrichment of causal variants within functional classes
  - Enables users to unbiasedly select from a (potentially) large pool functional annotations that are most phenotypically relevant
5. (optional) Model population-specific LD patterns.
6. (optional) Leverage cross
7. (optional) Improved computational performance via Importance Sampling. 

## Installation
The software has two dependencies: [1] Eigen v3.2 (matrix library) [2] NLopt v2.4.2 (optimization library) which are packaged with PAINTOR in order to simplify installation. Please see the [Eigen homepage](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [NLopt homepage](http://ab-initio.mit.edu/wiki/index.php/NLopt) for more information.

Download the latest [version](https://github.com/gkichaev/PAINTOR_FineMapping/releases) of the software into your desired target directory. Then unpack and install the software with the following commands:

`git clone https://github.com/gkichaev/PAINTOR_V3.0.git`

`cd PAINTOR_V3.0`

`bash install.sh`

This will create an executable "PAINTOR". Sample data is provided with the package. To test that the installation worked properly, type:

`./PAINTOR -input SampleData/input.files -in SampleData/ -out SampleData/ -Zhead Zscore -LDname ld -enumerate 2 -annotations DHS`

If everything worked correctly the final sum of log Bayes Factors should be: `653.892389`

For quick start simply type:

`./PAINTOR`

For detailed information on input files and command line flags see the user manual provided.

## Functional Annotations
We have compiled library of functional annotations that you may find useful. You can download it [here] (https://ucla.box.com/s/x47apvgv51au1rlmuat8m4zdjhcniv2d). This large compendium includes most of the Roadmap/ENCODE data as well as other regulatory and genic annotations. 
