# PAINTOR
Probabilistic Annotation INtegraTOR

## UPDATE 01/10/18
Announcing PAINTOR v3.1!
1. The new version updates the PAINTOR approximate inference scheme to use an efficient Gibbs sampling algorithm to sample directly from the posterior. Specify with the `-mcmc` flag. Note that exact inference (limited to k causals) can still be done with the `-enumerate [k]` flag.  
2. The prior effect size variance is estimated directly from the data rather than being fixed apriori. This allows us accomodate variability in effect sizes across fine-mapping regions. We use truncated SVD to estimate N*h2g for each locus. Use `-prop_ld` flag to modify the proportion of the LD spectrum to keep (default = 0.95).
3. Fixes bug in NLopt package that would result in "Optimization Errors" being thrown.
4. More logging of relevant output to aid debugging. 

## UPDATE 01/07/17
Announcing PAINTOR v3.0! The new version has enhancements that improve computational effiency, statistical robustness, as well as having expanded functionality to leverage multiple traits. In adddition, we have developed a visualiziation tool, [PAINTOR-CANVIS](https://github.com/gkichaev/PAINTOR_V3.0/tree/master/CANVIS), to produce publication-ready plots for the output of PAINTOR as seen below.

For legacy purposes, we leave available [PAINTOR 2.1](https://github.com/gkichaev/PAINTOR_V2.1), though we recommend using this latest version for most accurate results.

## Description

<img align="right" src="CANVIS/canvis.png">

We provide a command line implementation of the PAINTOR frameworks described in Kichaev et al. [(PLOS Genetics, 2014)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1004722),  [(American Journal of Human Genetics, 2015)](http://www.cell.com/ajhg/abstract/S0002-9297(15)00243-8), and  [(Bioinformatics, 2016)](http://bioinformatics.oxfordjournals.org/content/early/2016/10/16/bioinformatics.btw615).  Briefly, PAINTOR is a statistical fine-mapping method that integrates functional genomic data with association strength from potentially multiple populations (or traits) to prioritize variants for follow-up analysis. The software runs on multiple fine-mapping loci and/or populations/traits simultaneously and takes as input the following data for each set of SNPs at a locus


1. Summary Association Statistics (Z-scores)
2. Linkage Disequilibrium Matrix/Matrices (Pairwise Pearson correlations coefficients between each SNP)
3. Functional Annotation Matrix (Binary indicator of annotation membership (i.e. if entry {i,k} = 1, then SNP i is a member of annotation K).

#### Key Features

1. Outputs a probability for a SNP to be causal which can subsequently be used to prioritize variants
2. Can model multiple causal variants at any risk locus
3. Leverage functional genomic data as a prior probability to improve prioritization
  - This prior probability is not pre-specified, but rather, learned directly from the data via Empirical Bayes.
4. Quantify enrichment of causal variants within functional classes
  - Enables users to unbiasedly select from a (potentially) large pool functional annotations that are most phenotypically relevant
5. Fully Bayesian treatment of causal effect sizes
6. (optional) Model population-specific LD patterns when doing multi-ethnic fine-mapping.
7. (optional) Joint inference across traits when doing multi-trait fine-mapping.
8. (optional) Approximate inference via Gibbs Sampling.

#### For detailed information about input file formats, command line flags, and recommended analysis pipelines please see the [wiki](https://github.com/gkichaev/PAINTOR_V3.0/wiki)

## Installation
The software has two dependencies: [1] Eigen v3.2 (matrix library) [2] NLopt v2.4.2 (optimization library) which are packaged with PAINTOR in order to simplify installation. Please see the [Eigen homepage](http://eigen.tuxfamily.org/index.php?title=Main_Page) and [NLopt homepage](http://ab-initio.mit.edu/wiki/index.php/NLopt) for more information. Note that compiling requires gcc V4.9 (or greater).

For quick installation:

`git clone https://github.com/gkichaev/PAINTOR_V3.0.git`

`cd PAINTOR_V3.0`

`bash install.sh`

This will create an executable "PAINTOR". Sample data is provided with the package. To test that the installation worked properly, type:

`./PAINTOR -input SampleData/input.files -in SampleData/ -out SampleData/ -Zhead Zscore -LDname ld -enumerate 2 -annotations DHS`

If everything worked correctly the final sum of log Bayes Factors should be: `658.648`

For quick start simply type:

`./PAINTOR`


## Functional Annotations
We have compiled library of functional annotations that you may find useful. This large compendium includes .bed files for most of the Roadmap/ENCODE data as well as other regulatory and genic annotations.  Please see the  [wiki] (https://github.com/gkichaev/PAINTOR_V3.0/wiki/2b.-Overlapping-annotations) for more information and download link.
