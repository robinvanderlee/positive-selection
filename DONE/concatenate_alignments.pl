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
use Bio::Align::Utilities;

require('functions.pl');


### output directory for sequence files
my $OUTPUT_ROOT = "sequences";
my $OUTPUT_CDS = "cds";
my $OUTPUT_MASKED_ALIGNMENTS = "$OUTPUT_ROOT/prank-codon-masked";


### 
my $BASE_DIRECTORY = "$OUTPUT_MASKED_ALIGNMENTS";
opendir(BASEDIR, $BASE_DIRECTORY) or die $!;


###########################
####### MAIN SCRIPT #######
###########################
print STDERR "Concatenating the alignments in $BASE_DIRECTORY into one alignment...\n\n";
my $count = 0;
my @array_of_alignments = ();

### loop over all sorted alignments results file in the base directory
while (my $basedir_file = readdir(BASEDIR)) {
  next unless(-f "$BASE_DIRECTORY/$basedir_file"); # only loop over files
  
  # only loop over relevant files, i.e. the # GUIDANCE & TCS MASKED ALIGNMENT, SORTED BY SPECIES
  if($basedir_file =~ /^(ENSG\d+)__cds\.prank\-codon\-guidance\-tcs\-masked\-species\-sorted\.aln\.fa$/){

    ### PRINT COUNTER
    $count++;
    print STDERR $count . " MSAs read\r" if (($count % 100) == 0);

    my $current_ensembl_id = $1;

    ### GET ALIGNMENT
    my $current_alignment = undef;
    eval { $current_alignment = readFastaFile("$OUTPUT_MASKED_ALIGNMENTS" . "/" . $basedir_file) }; # GUIDANCE & TCS MASKED ALIGNMENT, SORTED BY SPECIES
    die if ($@);

    # print $basedir_file . "\t" . $current_ensembl_id . "\t" . $current_alignment->num_sequences . "\n";
    
    ### push to alignment array
    push(@array_of_alignments, $current_alignment);
  }

}
print STDERR $count . " MSAs read\n";
print scalar(@array_of_alignments) . " alignments in array\n";
closedir BASEDIR;



### concatenate alignments into one alignment
 # Title     : cat
 # Usage     : $aln123 = cat($aln1, $aln2, $aln3)
 # Function  : Concatenates alignment objects. Sequences are identified by id.
 #             An error will be thrown if the sequence ids are not unique in the
 #             first alignment. If any ids are not present or not unique in any
 #             of the additional alignments then those sequences are omitted from
 #             the concatenated alignment, and a warning is issued. An error will
 #             be thrown if any of the alignments are not flush, since
 #             concatenating such alignments is unlikely to make biological
 #             sense.
 # Returns   : A new Bio::SimpleAlign object
 # Args      : A list of Bio::SimpleAlign objects
#my $concatenated_alignment = Bio::Align::Utilities->cat($array_of_alignments[1], $array_of_alignments[10]); ### ERROR Can't locate object method "id" via package "Bio::Align::Utilities" at /mnt/datas/rvdlee/lib/lib/perl5/Bio/Align/Utilities.pm line 372.

my %sequences = ();

foreach my $alignment ( @array_of_alignments ) {
  foreach my $seq ($alignment->each_seq) {
    my $identifier = $seq->display_id;
    my $species = "";
    $species = $1 if($identifier =~ /^(\d+__\S+?)__/);

    $sequences{$species} .= $seq->seq; # append the sequence
  }
}


### write sequences
my $outfile = "$OUTPUT_ROOT/concatenated_alignment__9primates__cds.prank-codon-guidance-tcs-masked.aln.fa";
my $seqout = Bio::SeqIO->new
                            (-file => ">$outfile",
                             -format => 'fasta');

while( my ($seqname,$seqstr) =  each %sequences )  {
  my $seq = Bio::Seq->new(-id => $seqname, -seq => $seqstr);
  $seqout->write_seq($seq);
}

print "File $outfile has been generated\n";



### FINISH
print STDERR "\nFINISHED\n";






#########################
####### FUNCTIONS #######
#########################
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
