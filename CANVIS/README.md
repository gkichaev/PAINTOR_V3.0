# CANVIS (Correlation ANnotation Visualization)

## Description
We provide here documentation for CANVIS, a visualization tool utilized for the output of PAINTOR. We will describe the  neccessary libraries, file formats, and input parameters, as well as a sample script. 

The final visualization is composed of the following:

1. Scatterplot of location versus posterior probabilites with a specified credible set
2. Annotation bars
3. Scatterplot of location versus -log10(pvalue)   
4. Correlation heatmap of LD matrix 

The following subplots are output into one PDF file designated by the  desired output parameter. 

## Installation
Download the latest software into target directory. The neccessary libraries can be installed on either OSX or Linux with `install_mac.sh` or `install_linux.sh`. Additionally, the madatory dependencies are as follows:

- Python 2.7
- numpy
- scipy
- matplotlib

Sample data is also provided in a folder, `CANVIS_Sample`. To test the installation, run the following command inside of the sample folder:

```
$ python CANVIS.py -l chr4.3473139.rs6831256.post.filt.300 -z ldl.Zscore tg.Zscore -r chr4.3473139.rs6831256.ld.filt.300 -a chr4.3473139.rs6831256.annot.filt.300 -s E066.H3K27ac.narrowPeak.Adult_Liver E066.H3K4me1.narrowPeak.Adult_Liver -t 99 -i 3381705 3507346

```

## File Formats
The input files for CANVIS are the output files from PAINTOR, and these are described in more detail below. 

#### Locus File
The locus file is space delimted, and the first row is expected to be a header row containing the labels of each column. This file must contain at least one column of zscores with the name of the specific zscore as the heading, one column of posterior probabilites with the heading `Posterior_Prob`, and a column holding the positions with the heading `pos`. CANVAS can plot between 1 - 3 zscores. An example of a correctly formatted file is shown below: 

```
chr pos rsid hdl.A0 hdl.A1 hdl.Zscore ldl.Zscore tc.Zscore tg.Zscore Posterior_Prob
chr4 3349626 rs74823150 A G -0.594734 0.696424 0.507638 0.345514 2.7728e-06
chr4 3350248 rs2749779 A G -1.320755 1.068966 1.754386 2.153846 1.87298e-06

```
#### Annotation File
The annotation file is space delimited, and should contain a column per annotation where each position has a 0 or 1 to represent a specific annotation. The first row is assumed to be a header row containing the names of the annotations. The file can contain multiple columns, but a max of 5 annotations can be plotted. 

```
E066-H3K27ac.narrowPeak.Adult_Liver E066-H3K4me1.narrowPeak.Adult_Liver
0 1
0 0
1 0
```
#### LD file
The ld matrix file is a space delimited file with no header row or indexing column. Each row should contain the correlations corresponding to other locations of the sample. 

```
0.20794 -0.33251 -0.0018495 -0.33251 0.7164 -0.33181 0.70725 0.62996 -0.33181 -0.33181 -0.33181
```

## Command Line Flags  
Here are the following input parameters specified at the command line: 

- `--locus [-l]` path to file with fine-mapping locus
	- Mandatory file; if not included, program will exit
- `--zscores [-z]` specific zscores to be plotted (between 1-3)
	- Mandatory file; if not included, program will exit
- `--annotations [-a]` path the file with annotations 
	- Optional; if not included, no annotation bars will be included
- `--specific_annotations[-s]` specific annotations to be plotted (between 1-5)
	- list specific annotations to plot (recommended 3 max)
	- if not specified, program will use all annotations listed in the heading of the file
- `--ld_name [-r]` path to file with ld matrix
	- Optional; if not included, no ld matrix will be plotted
	- recommended < 300 entries; if there are more entries, program will display a warning message that the program might perform significantly slower
- `--threshold [-t]` threshold for credible set, an number (0,100); default: `0`
	- Optional
	- if user specifies an invalid threshold, then default is enabled
- `--greyscale [-g]` `y` or `n`, flag for figure in greyscale; default: `n`
- `--output [-o]` desired name of output file [default: fig_final]
- `--interval [-i]` designated interval [default: all locations]
	- if the right and/or left interval are out of bounds or if there are no SNPs in the specified interval, they will be replaced with the minimum or maximum location respectively






