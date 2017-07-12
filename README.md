# Genome-scale detection of positive selection in 9 primates

This page contains scripts, analyses and procedures developed in:

> **Genome-scale detection of positive selection in 9 primates predicts human-virus evolutionary conflicts**
> Robin van der Lee, Laurens Wiel, Teunis J.P. van Dam, Martijn A. Huynen  
> bioRxiv 131680; https://doi.org/10.1101/131680

> Centre for Molecular and Biomolecular Informatics - Radboud Institute for Molecular Life Sciences<br/>
> Radboud university medical center<br/>
> Nijmegen, The Netherlands

Please consider citing the paper if you found this resource useful.


## Issues & Contact

Note that these scripts were not designed to function as a fully-automated pipeline, but rather as a series of steps with extensive manual quality control between them. It will therefore not be straightforward to run all steps smoothly in one go. Feel free to contact me if you run into any issue with individual steps.

> **robinvanderlee AT gmail DOT com**<br/><br/>
> [Google Scholar](https://scholar.google.co.uk/citations?user=ISYCcUUAAAAJ)<br/>
> [ORCID](http://orcid.org/0000-0001-7391-9438)<br/>
> [Twitter](https://twitter.com/robinvdlee)<br/>
> [LinkedIn](http://nl.linkedin.com/in/robinvdlee)


## Requirements

Scripts depend on various programs and modules to run.

###### Perl
- Install Perl 5: https://www.perl.org/
- Install various modules using `cpan` or `cpanm`
	- `cpanm DBI`
	- `cpanm DBD::mysql`  (requires a working mysql installation, https://www.mysql.com/)
- Download the various helper scripts that are also part of this GitHub repository
	- `functions.pl`
	- Scripts in `Ensembl_API_subroutines`

###### BioPerl
- Install BioPerl `cpanm Bio::Perl` (http://bioperl.org/)

###### Ensembl API
- Install the Ensembl API (http://www.ensembl.org/info/docs/api/index.html)

###### R
- Install R: https://www.r-project.org/
- Install the following packages:
	- ``

R
	+ modules
	?? gebruik ik dit script? ./Biological_correlations_PSR/basic_statistics_of_data/aa_freqs/analyze_aa_freqs.r:require(lattice)
library("Biostrings")
library("biomaRt")
library("ggplot2")
library("phangorn")
library("session")
library(GenomicRanges)
library(MASS)
library(data.table)
library(dplyr)
library(ggplot2)
library(gplots)
library(gridExtra)
library(jsonlite)
library(parallel)
library(plotrix)
library(rtracklayer)
library(scales)
library(vioplot)





Parallel
PRANK
PAML/CODEML
	make sure they are in your path

Jalview
...

pal2nal?
muscle
mafft

GUIDANCE
t_coffee

NOTE THAT ENSEMBL VERSION USED FOR THE PAPER IS XXX
Ensembl release 78, December 2014 (http://dec2014.archive.ensembl.org/)


## Steps

Please see the `Materials and Methods` section of the paper for theory and detailed explanations.


### 1. One-to-one orthologs
Obtain one-to-one ortholog clusters for nine primates with high-coverage whole-genome sequences. These scripts can be edited to obtain orthology clusters for (i) a different set of species than the ones we use here, and (ii) different homology relationships than the one-to-one filter we use.<br/>
Two methods, same result:

###### 1a. Ensembl API method
1. Fetch orthology information using the Ensembl API: `get_one2one_orthologs__Ensembl_API.pl`
2. Clean the results using: `get_one2one_orthologs__Ensembl_API__process_orthologs.r`

###### 1b. Ensembl BioMart method 
1. From BioMart, first get orthology information for all species of interest
![alt text](Images/Step1b__1.png)
![alt text](Images/Step1b__2.png)
2. Combine the acquired ortholog information using `get_one2one_orthologs__combine_biomart_orthology_lists.r`


### 2. Sequences
For all one-to-one orthologs, get the coding DNA (cds) and the corresponding protein sequences from the Ensembl Compara gene trees (as these trees are the basis for the orthology calls).

###### 2a. Get sequences
1. Run `start_parallel_sequence_fetching.pl`, which reads the ortholog cluster information, divides the genes in batches of 1000 (`$batchsize` can be changed in the script) and prints the instructions for the next step:
2. `get_sequences_for_one2one_orthologs_from_gene_tree_pipeline.pl`. This step fetches the sequences. E.g.
```
STARTING PARALLEL INSTANCE 11
from 10001 to 11000
1000 genes
screen -S i11
perl get_sequences_for_one2one_orthologs_from_gene_tree_pipeline.pl -p -f sequences/parallel_instance_11.txt
```
Note that these steps also fetch the alignments underlying the Compara gene trees, filtered for the species of interest. These are not required for later steps.

###### 2b. Check sequences
`check_compatibility_protein_cds_sequences.pl`. Tends to only complain at annotated selenocysteine residues.


### 3. Alignments


###### 3a. PRANK codon-based multiple alignment
We  codon-based alignments nucleotide alignments

. The
PRANK algorithm to some extent prevents alignment of nonhomologous regions by flagging
gaps made during different stages of progressive alignment and permitting their reuse without
further penalties (Löytynoja and Goldman 2008). As PRANK has been shown to provide
better initial alignments for large-scale positive selection detection (Schneider et al. 2009;
Fletcher and Yang 2010; Markova-Raina and Petrov 2011; Jordan and Goldman 2012;
Privman et al. 2012; Moretti et al. 2014), we obtained multiple alignments of the primate
ortholog clusters using the PRANK codon mode (prank +F –codon; v.140603). 
We used the
default settings of (i) obtaining a guide tree from MAFFT for the progressive alignment
procedure and (ii) selecting the best alignment from five iterations


Some useful command to monitor progress:
```
XX
```



###### 3b. GUIDANCE - assessment and masking



###### 3c. TCS - assessment and masking



###### 3d. Translate alignments
Steps are all CDS/codon based, but this is handy for checking alignments, analyzing results etc






## Supplementary data and material

Additional data and material can be found at:
- This GitHub page: [Supplementary data and material](Supplementary_data_and_material)<br/>

and
- http://www.cmbi.umcn.nl/~rvdlee/positive_selection/






