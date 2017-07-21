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
use Bio::AlignIO;


### command line arguments
my $phylip_alignment_file = "";
my $template_ctl_file = "";
if(scalar @ARGV != 2){
	die "Please provide a phylip file (argument 1), and a template codeml ctl file (argument 2)\n";
} else {
	$phylip_alignment_file = shift @ARGV;
	die "Please provide a phylip file\n" if $phylip_alignment_file !~ /\.aln\.phy$/;
	
	$template_ctl_file = shift @ARGV;
	die "Please provide a template codeml ctl file\n" if $template_ctl_file !~ /template\.ctl$/;
}


### extract ensembl identifier
my $ensembl_identifier = "";
if($phylip_alignment_file =~ /(ENSG\d+)__/){
	# print $1 . "\n";
	$ensembl_identifier = $1;
} else {
	die "Cannot extract ENSG\\d+ identifier from the file name: $phylip_alignment_file\n";
}


### create directory structure
my $BASEDIR = "codeml_results";
`mkdir $BASEDIR` unless -d $BASEDIR; # unless it already exists
my $NEW_WORKING_DIR = "$BASEDIR/$ensembl_identifier";
`mkdir $NEW_WORKING_DIR` unless -d $NEW_WORKING_DIR; # unless it already exists

### change template .ctl to deal with the provided alignment file

# new ctl file
my $specific_ctl_file = "";
if ($template_ctl_file =~ /^(\S+?)\_\_/){
	$specific_ctl_file = $ensembl_identifier . "_" . $1 . ".ctl";
} else {
	die "Cannot create specific ctl file from: $template_ctl_file\n";
}
# print $specific_ctl_file . "\n";

open(TEMPLATE_CTL_IN, "<$template_ctl_file");
open(SPECIFIC_CTL_OUT, ">$NEW_WORKING_DIR/$specific_ctl_file");

# read and replace
while(my $line = <TEMPLATE_CTL_IN>){
	if($line =~ /<ALIGNMENT_FILE_REPLACE>/){
		$line =~ s/<ALIGNMENT_FILE_REPLACE>/\.\.\/\.\.\/$phylip_alignment_file/;
	}

	if($line =~ /<OUTFILE_REPLACE>/){
		$line =~ s/<OUTFILE_REPLACE>/$ensembl_identifier/;
	}
	print SPECIFIC_CTL_OUT $line;
}

close TEMPLATE_CTL_IN;
close SPECIFIC_CTL_OUT;


### change working directory and start up codeml
chdir $NEW_WORKING_DIR;
`codeml $specific_ctl_file > $ensembl_identifier.screen_output`;

if($? == -1){
	die "codeml failed: $!\n";
} else {
	exit 0;
}
