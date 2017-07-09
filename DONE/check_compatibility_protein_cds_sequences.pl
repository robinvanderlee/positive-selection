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
use Bio::Seq;
use Bio::SeqIO;

### output directory for sequence files
my $OUTPUT_ROOT = "sequences";
my $OUTPUT_PROTEIN = "protein";
my $OUTPUT_CDS = "cds";

### 
my $directory = "$OUTPUT_ROOT/$OUTPUT_CDS";
opendir(DIR, $directory) or die $!;

### which sequences to compare?
my $file_suffix = "fa";
if(defined($ARGV[0])){
	if($ARGV[0] eq "-compara"){
		$file_suffix = "compara.aln.fa";
	}
}




####### SCRIPT
print "Checking protein and cds sequence compatibility...\n\n";

my $count = 0;
my $count_total_seq = 0;
my $problem_count = 0;
my %problem_ortholog_groups = ();
while (my $file = readdir(DIR)) {
	next unless $file =~ /^(.+?)__cds.$file_suffix$/;
	my $id = $1;

	my $proteinfasta_io  = Bio::SeqIO->new(-file => "<$OUTPUT_ROOT/$OUTPUT_PROTEIN/$id\_\_protein.$file_suffix",
	                            	   -format => 'fasta');
	my @protein_seqs = ();
	while ( my $seq = $proteinfasta_io->next_seq() ) {
		push(@protein_seqs, $seq);
	}


	my $cdsfasta_io  = Bio::SeqIO->new(-file => "<$directory/$file",
	                           		   -format => 'fasta');

	while ( my $seq = $cdsfasta_io->next_seq() ) {
	  my $translated_cds_seq = $seq->translate(-terminator => '*', -throw => 1)->seq;
	  # my $translated_cds_seq = $seq->translate(-complete => 1, -throw => 0)->seq;
	  $translated_cds_seq =~ s/\*$//; # remove terminator

	  foreach my $protein_seq (@protein_seqs){
	  	if($protein_seq->id eq $seq->id){
	  		if($protein_seq->seq ne $translated_cds_seq){
	  			print "Translated cds sequence is different from protein sequence:\n" . $id . "\t" . $seq->id . "\n";
	  			print ">Protein sequence\n" . $protein_seq->seq . "\n";
	  			print ">Translated cds sequence\n" . $translated_cds_seq . "\n";
	  			print ">cds sequence\n" . $seq->seq . "\n";
	  			print "\n==\n\n";
	  			$problem_count++;
	  			$problem_ortholog_groups{$id} = 1;
		  	}
		  	last;
	  	}
	  }

	  
	  $count_total_seq++;
	}

	$count++;
	print STDERR $count . " sets of orthologous sequences done\r" if (($count % 100) == 0);
}
print STDERR $count . " sets of orthologous sequences done\r";

print "\nAnalyzed $count sets of orthologous sequences, containing $count_total_seq sequences\n";
print scalar(keys %problem_ortholog_groups) . " problem cases, with $problem_count problematic sequences\n";
print "FINISHED\n";
