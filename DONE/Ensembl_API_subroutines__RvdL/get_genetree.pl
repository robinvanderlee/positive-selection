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

### Fetches the GeneTree object for an Ensembl Genemember
sub get_genetree{
	# ARGS:
	# GeneTree adaptor
	my $genetree_adaptor = shift @_;
	
	# GeneMember object
	my $member = shift @_;

	# fetches from the database all the gene trees that contains this member
	my $genetree = $genetree_adaptor->fetch_all_by_Member($member)->[0];
	die "Cannot fetch GeneTree for " . $member->stable_id . "\n" unless (defined $genetree);

	return $genetree;
}

return 1;
