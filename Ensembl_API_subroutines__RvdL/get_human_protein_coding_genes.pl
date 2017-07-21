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

sub get_human_protein_coding_genes{
	# ARGS:
	# 1. Ensembl API registry

	my $reg = $_[0];
	
	# Get adaptors
	my $gene_adaptor = $reg->get_adaptor( 'Human', 'Core', 'Gene' );

	print "Getting all human protein coding genes...\n";

	# fetch all human genes
	my @genes = @{ $gene_adaptor->fetch_all() };

	# filter for protein coding genes
	my @protein_coding_genes = ();
	foreach my $gene (@genes){
		if($gene->biotype =~ /protein_coding/){
			push(@protein_coding_genes, $gene->stable_id);
		}
	}

	return @protein_coding_genes;
}

return 1;
