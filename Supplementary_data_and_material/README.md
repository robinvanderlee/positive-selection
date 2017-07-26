# Supplementary Data and Material

This directory contains Supplementary Data and Material for

> **Genome-scale detection of positive selection in 9 primates predicts human-virus evolutionary conflicts**
> Robin van der Lee, Laurens Wiel, Teunis J.P. van Dam, Martijn A. Huynen  
> bioRxiv 131680; https://doi.org/10.1101/131680

> Centre for Molecular and Biomolecular Informatics - Radboud Institute for Molecular Life Sciences<br/>
> Radboud university medical center<br/>
> Nijmegen, The Netherlands

> Correspondence: robinvanderlee AT gmail DOT com

If you found this resource useful, please consider citing the paper.


## Original coding and proteins sequences used in our analyses
- Protein-coding DNA sequences<br/>
`cds_sequences_11170_genes.tar.gz`<br/>
(11,170 files, 9 sequences per file, 1 per primate)

- Protein sequences<br/>
`protein_sequences_11170_genes.tar.gz`<br/>
(11,170 files, 9 sequences per file, 1 per primate)


## PRANK alignments, filtered and masked using GUIDANCE and TCS
- Protein-coding DNA sequence alignments<br/>
`prank-codon-guidance-tcs-masked-species-sorted.aln__11096.tar.gz`<br/>
`N`: undetermined nucleotides<br/>
`n`: masked nucleotides

- Translated alignments<br/>
`prank-codon-guidance-tcs-masked-species-sorted.aln.translated__11096.tar.gz`<br/>
`X`: undetermined amino acids<br/>
`o`: masked amino acids

- Mapping tables of the alignment column to the position in the human sequence<br/>
`prank-codon-guidance-tcs-masked-species-sorted.aln.translated.mapping-table__11096.tar.gz`<br/>
Column 1: human amino acid<br/>
Column 2: columns number of the alignment<br/>
Column 3: human sequence position<br/>

- Concatenation of the 11,096 protein-coding DNA sequence alignments<br/>
`concatenated_alignment__11096_genes__9primates__cds.prank-codon-guidance-tcs-masked-species-sorted.aln.fa.gz`
`concatenated_alignment__11096_genes__9primates__cds.prank-codon-guidance-tcs-masked-species-sorted.aln.phy.gz`


## Visualization and quality control of alignments and positive selection profiles
- Barcode plots of positive selection profiles<br/>
  * Just the curated set of 331 positively selected genes (PSG): `barcode_plots__331_PSG.tar.gz`<br/>
  * All 416 original, apparent positive selected genes (aPSG): `barcode_plots__416_aPSG.tar.gz`<br/>

- Jalview alignment annotations<br/>
`jalview_alignment_annotations__331_PSG.tar.gz`<br/>
`jalview_alignment_annotations__416_aPSG.tar.gz`

- Jalview sequence annotations<br/>
`jalview_sequence_annotations__331_PSG.tar.gz`<br/>
`jalview_sequence_annotations__416_aPSG.tar.gz`<br/>

<br/>

Jalview annotation files contain annotations for the positively selected residues (PSR), as well as exon coordinates mapped to protein sequences. They can be loaded together with the correct corresponding alignment using (on Mac):
```
Java	-Djava.ext.dirs=/Applications/Jalview/lib/
	-cp /Applications/Jalview/jalview.jar jalview.b.Jalview
	-open <ENSEMBL_ID>__cds.prank-codon-guidance-tcs-masked-species-sorted.aln.translated.fa
	-annotations <ENSEMBL_ID>.jalview_aln_feature
	-features <ENSEMBL_ID>.jalview_seq_feature
```

## Configuration files for PAML codeml
- Inference of reference phylogenetic trees<br/>
`codeml_M0_F3X4_tree.ctl`<br/>
`codeml_M0_F61_tree.ctl`<br/>

- Inference of positive selection<br/>
`codeml_M1avM2a_F3X4__large-scale-analysis__template.ctl`<br/>
`codeml_M1avM2a_F61__large-scale-analysis__template.ctl`<br/>
`codeml_M7vM8_F3X4__large-scale-analysis__template.ctl`<br/>
`codeml_M7vM8_F61__large-scale-analysis__template.ctl`<br/>


## Phylogenetic trees
- codeml reference phylogenetic trees<br/>
`codeml_M0_tree__unrooted_tree__F3X4.tre`<br/>
`codeml_M0_tree__unrooted_tree__F61.tre`<br/>

- Other phylogenetic trees<br/>
`Ensembl78__9primates__with_taxon_id__unrooted.tre`<br/>
`Perelman_et_al__9primates__unrooted.tre`<br/>
`RAxML_bestTree__9primates__unrooted.tre`<br/>


## GUIDANCE source code fix for running PRANK
- `guidance_prank_fix.pdf`
- `Guidance__prank_fix.pm`
