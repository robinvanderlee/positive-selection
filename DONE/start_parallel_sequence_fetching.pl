#!/usr/bin/env perl
#####################################
## Robin van der Lee               ##
## robinvanderlee AT gmail DOT com ##
############################################################################################################
## Genome-scale detection of positive selection in 9 primates predicts human-virus evolutionary conflicts ##
## Robin van der Lee, Laurens Wiel, Teunis J.P. van Dam, Martijn A. Huynen                                ##
############################################################################################################
use warnings;
use strict;
use File::Path qw( make_path );

my $batchsize = 1000;

### output directory
my $OUTPUT_ROOT = "sequences";
my $OUTPUT_PROTEIN = "protein";
my $OUTPUT_CDS = "cds";

my @dirs = ("$OUTPUT_ROOT/$OUTPUT_PROTEIN", "$OUTPUT_ROOT/$OUTPUT_CDS");
foreach my $dir(@dirs){
  if(!-d $dir){
    make_path($dir) or die "Cannot create output directories";
  }
}

###################################################################################################
### read one2one orthologs
# Ensembl.Gene.ID Chimpanzee.Ensembl.Gene.ID  Gibbon.Ensembl.Gene.ID  Gorilla.Ensembl.Gene.ID Marmoset.Ensembl.Gene.ID  Olive.baboon.Ensembl.Gene.ID  Orangutan.Ensembl.Gene.ID Vervet.AGM.Ensembl.Gene.ID  Macaque.Ensembl.Gene.ID
# ENSG00000000003 ENSPTRG00000022087  ENSNLEG00000007505  ENSGGOG00000010003  ENSCJAG00000004288  ENSPANG00000018367  ENSPPYG00000020538  ENSCSAG00000010385  ENSMMUG00000008209
# ENSG00000000005 ENSPTRG00000022086  ENSNLEG00000007504  ENSGGOG00000017007  ENSCJAG00000004270  ENSPANG00000016641  ENSPPYG00000020537  ENSCSAG00000010392  ENSMMUG00000008212
# ENSG00000000419 ENSPTRG00000013625  ENSNLEG00000006998  ENSGGOG00000008993  ENSCJAG00000000551  ENSPANG00000000450  ENSPPYG00000011131  ENSCSAG00000015216  ENSMMUG00000002759
# ...
my $ortholog_file = "biomart_one2one_ortholog_clusters.txt";
my @total_genes_list = ();
open(IN, "<$ortholog_file") or die "Cannot open $ortholog_file\n";
while(<IN>){
	next if $. == 1;
	chomp;
	my @F = split /\t/;
	push(@total_genes_list, $F[0]);
};
my $total_genes = scalar(@total_genes_list); # minus header
close IN;



###################################################################################################
my $instances = int($total_genes / $batchsize) + 1;

for(my $i = 1 ; $i <= $instances; $i++){

	my $from = 1 + ( ($i - 1) * $batchsize);
	my $to = $from + $batchsize - 1;
	$to = $total_genes if $to > $total_genes;

	print "STARTING PARALLEL INSTANCE $i\nfrom $from to $to\n";

	my @total_genes_list_selected = @total_genes_list[($from - 1) .. ($to - 1)];
	my $total_genes_selected = scalar @total_genes_list_selected;

	print $total_genes_selected . " genes\n";

	my $outfile = $OUTPUT_ROOT . "/parallel_instance_" . $i . ".txt";
	my $outputfile = $OUTPUT_ROOT . "/parallel_instance_" . $i . ".txt_output";
	open(OUT, ">$outfile");
	foreach(@total_genes_list_selected){
		print OUT $_ . "\n";
	}
	close OUT;

	my $cmd2 = "screen -S i$i";
    print "$cmd2\n";

	#my $cmd = "perl get_sequences_for_one2one_orthologs_from_gene_tree_pipeline.pl -p -f $outfile > $outputfile &";
	my $cmd = "perl get_sequences_for_one2one_orthologs_from_gene_tree_pipeline.pl -p -f $outfile";
	print "$cmd\n";

	# system $cmd;

	print "\n";
}


