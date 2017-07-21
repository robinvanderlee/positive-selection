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
use Bio::AlignIO;

require('functions.pl');

### command line arguments
# 1: infile
my $fasta_alignment_file = "";
if(scalar @ARGV == 0){
	die "Please provide a fasta input file\n";
} else {
	$fasta_alignment_file = shift @ARGV;
	die "Please provide a fasta input file\n" if $fasta_alignment_file !~ /fa[s]?$/;
}

# 2: outfile
my $outfile = "";
if(scalar @ARGV == 0){
	die "Please provide a fasta output file\n";
} else {
	$outfile = shift @ARGV;
}


### read fasta file
my $alignment_object = undef;
eval { $alignment_object = readFastaFile($fasta_alignment_file) };
die if ($@);

### define species of interest
my %species_of_interest = (
                        9606  =>  "homo_sapiens",
                        9598  =>  "pan_troglodytes",
                        9595  =>  "gorilla_gorilla",
                        9601  =>  "pongo_abelii",
                        61853 =>  "nomascus_leucogenys",
                        9544  =>  "macaca_mulatta",
                        9555  =>  "papio_anubis",
                        60711 =>  "chlorocebus_sabaeus",
                        9483  =>  "callithrix_jacchus"
                      ); # 9 Simian primates
my @species_of_interest_display_order = (9606, 9598, 9595, 9601, 61853, 9544, 9555, 60711, 9483);

### sort sequences within alignments by species, based on taxonomy
my @sequence_objects_alignment_ordered = order_sequences_by_taxon($alignment_object, \@species_of_interest_display_order);

### write alignment to file
# my $outfile = $fasta_alignment_file;
# $outfile =~ s/(\.aln\.fa)$/-species-sorted$1/;
my $seqout = Bio::SeqIO->new
                            (-file => ">$outfile",
                             -format => 'fasta');

$seqout->write_seq(@sequence_objects_alignment_ordered);
# print "File $outfile has been generated\n";
# print STDERR "\nFINISHED\n";




#########################
####### FUNCTIONS #######
#########################
###########################################################################
###########################################################################
# sub writeAlignment{
# 	my ($alignment_file, $alignment) = @_;

# 	my $AlignIO_out = Bio::AlignIO->new(-file   => ">$alignment_file", #-fh     => \*STDOUT ,
# 	 									-format => 'fasta');
# 	$AlignIO_out->write_aln($alignment);
# 	# $AlignIO_out->write_aln($alignment->slice(1,1000,1));
# }

sub order_sequences_by_taxon{
  # first argument: alignment object
  my $alignment_object = shift @_;

  # second argument: array of taxonomy numbers
  my $species_order_ref = shift @_;
  my @species_order = @{$species_order_ref};

  # printArray(@species_order);

  # get array of Bio::Seq objects, with -accession_number set as taxon
  my @sequence_objects_alignment = get_sequence_object_array_from_alignment_object($alignment_object);
  # foreach(@sequence_objects_alignment){
  # 	print $_->display_name . "\n";
  # 	print $_->accession_number . "\n\n";
  # }

  my @seq_objects_ordered = ();
  foreach my $species_order_taxon (@species_order){
    foreach my $seq_object (@sequence_objects_alignment){
      if($species_order_taxon == $seq_object->accession_number){ # taxon number is stored in accession_number!
        push(@seq_objects_ordered, $seq_object);
      }
    }
  }

  return @seq_objects_ordered;
}

##########################################################################
sub get_sequence_object_array_from_alignment_object{
	# first argument: alignment object
  	my $alignment_object = shift @_;

	my @sequence_objects_alignment = ();
	foreach my $alignment_sequence ($alignment_object->each_seq) {
	  my @F = split /\_/, $alignment_sequence->display_id;
	  my $alignment_sequence_taxon = $F[0];

	  my $sequence_object = Bio::Seq->new(
	                                -seq        => $alignment_sequence->seq,
	                                -display_id => $alignment_sequence->display_id,
	                                -accession_number  => $alignment_sequence_taxon, # store taxon in accession_number
	                            );
	  push(@sequence_objects_alignment, $sequence_object)
	}

	# returns:  array of Bio::Seq objects, with -accession_number set as taxon
	return @sequence_objects_alignment;
}

###########################################################################
sub readFastaFile{
	my $file = shift @_;
	
	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#read alignment file
	my $AlignIO_in  = Bio::AlignIO->new(-file   => $file ,
                                	    -format => 'fasta');
    my $alignment = $AlignIO_in->next_aln(); # Returns : a Bio::Align::AlignI compliant object
	$alignment->set_displayname_flat(); # Function  : Makes all the sequences be displayed as just their name, not name/start-end

    #throw error if file contains another alignment
    my $alignment2 = $AlignIO_in->next_aln();
    if(defined $alignment2){
      die "$file contains multiple alignments";
    }

    return $alignment;
}

###########################################################################
sub checkIfExists{
	my $file = shift @_;
	if(! -f $file){
		die "File $file does not exist";
	}
}
