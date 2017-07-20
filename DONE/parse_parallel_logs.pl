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

my $logfile = "all";
if(scalar @ARGV == 0){
	# die "Please provide a parallel log file, or specify \'all\'\n";
} else {
	$logfile = shift @ARGV;
}

if($logfile eq "all"){
	opendir(BASEDIR, ".") or die $!;

	while (my $file = readdir(BASEDIR)) {
		if($file =~ /\.log$/){
			print $file . "\n\n";
			parse_log_file($file);
			print "=" x 50 . "\n";
		}
	}

	closedir BASEDIR;
} else {
	parse_log_file($logfile);
}



########
sub parse_log_file {
	my $logfile = shift @_;

	open(IN, "<$logfile");

	# Sat Apr 18 14:39:08:rvdlee@fm:~/analysis_positive_selection_primates__FM$ head parallel_reverse-translate-prank-codon-translated-alignments.log 
	# Seq	Host	Starttime	JobRuntime	Send	Receive	Exitval	Signal	Command
	# 1	:	1429359852.252	     0.559	0	0	0	0	echo TEST
	# 2	:	1429359852.258	     0.564	0	0	0	0	echo TEST
	my $header = "";
	while(my $line = <IN>){
		chomp $line;
		
		if($. == 1){
			$header = $line;
			next;
		}

		my @F = split /\s+/, $line;

		my $job_no = shift @F;
		shift @F;
		shift @F;
		shift @F;
		shift @F;
		shift @F;
		my $exitval = shift @F;
		my $signal = shift @F;
		my $cmd = join(" ", @F);

		print $line . "\n" if $exitval != 0;
	}
	close IN;
}
