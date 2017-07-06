# Genome-scale detection of positive selection in 9 primates
TEST
This page contains scripts, analyses and procedures developed in:

> **Genome-scale detection of positive selection in 9 primates predicts human-virus evolutionary conflicts**
> Robin van der Lee, Laurens Wiel, Teunis J.P. van Dam, Martijn A. Huynen  
> bioRxiv 131680; https://doi.org/10.1101/131680

> Centre for Molecular and Biomolecular Informatics - Radboud Institute for Molecular Life Sciences<br/>
> Radboud university medical center<br/>
> Nijmegen, The Netherlands

If you found this resource useful, please consider citing the paper.


## Issues & Contact

Please note that these scripts were not designed to function as a fully-automated pipeline, but rather as a series of steps with extensive manual quality control between them. It may therefore not be straightforward to run all steps smoothly in one go. Feel free to contact me if you run into any issue with individual steps.

> Correspondence: robinvanderlee AT gmail DOT com

[] contact, scholar, researchgate, twitter, ....................


## Requirements

Scripts depend on various programs and modules to run.

###### Perl
- Install Perl 5: https://www.perl.org/
- Install various modules using `cpan` or `cpanm`
	- `cpanm DBI`
	- `cpanm DBD::mysql`  (requires a working mysql installation, https://www.mysql.com/)

###### BioPerl
- Install BioPerl   `cpanm Bio::Perl` (http://bioperl.org/)

###### Ensembl API
- Install the Ensembl API (http://www.ensembl.org/info/docs/api/index.html)




in PERL5PATH
		functions.pl 
		+ 	Ensembl_API_subroutines__RvdL
require('../Ensembl_API_subroutines__RvdL/get_genetree.pl');
require('../Ensembl_API_subroutines__RvdL/get_human_protein_coding_genes.pl');
require('Ensembl_API_subroutines__RvdL/get_genetree.pl');
require('Ensembl_API_subroutines__RvdL/get_human_protein_coding_genes.pl');
require('functions.pl');




R
	+ modules



Parallel
PRANK
PAML/CODEML
	make sure they are in your path

Jalview
...


NOTE THAT ENSEMBL VERSION USED FOR THE PAPER IS XXX
ENSEMBL 78



## Steps

Please see the `Materials and Methods` section of the paper for detailed explanations.

#### 1. One-to-one orthologs

Obtain one-to-one ortholog clusters for nine primates with high-coverage whole-genome sequences. These scripts can be edited to obtain orthology clusters for (i) a different set of species than the ones we use here, and (ii) different homology relationships than the one-to-one filter we use.**
Two methods, same result:

###### 1a. Ensembl API method
1. Fetch orthology information using the Ensembl API: `get_one2one_orthologs__Ensembl_API.pl`
2. Clean the results using: `get_one2one_orthologs__Ensembl_API__process_orthologs.r`

###### 1b. Ensembl BioMart method 
1. From BioMart, first get orthology information for all species of interest
![alt text](Images/Step1b__1.png)
![alt text](Images/Step1b__2.png)
2. Combine the acquired ortholog information using `get_one2one_orthologs__combine_biomart_orthology_lists.r`





## Supplementary data and material

Additional data and material can be found at:
- This GitHub page: [Supplementary data and material](Supplementary_data_and_material)**
and
- http://www.cmbi.umcn.nl/~rvdlee/positive_selection/






