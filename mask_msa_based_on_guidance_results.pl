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
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use File::Path qw( make_path );

require('functions.pl');





### output directory for sequence files
my $OUTPUT_ROOT = "sequences";
my $OUTPUT_CDS = "cds";
my $OUTPUT_GUIDANCE_PRANK_CODON = "$OUTPUT_ROOT/guidance-prank-codon";
my $OUTPUT_GUIDANCE_PRANK_CODON_ERRORS = "$OUTPUT_ROOT/guidance-prank-codon-errors";
if(!-d $OUTPUT_GUIDANCE_PRANK_CODON_ERRORS){
	make_path($OUTPUT_GUIDANCE_PRANK_CODON_ERRORS) or die "Cannot create output directory";
}
my $OUTPUT_MASKED_ALIGNMENTS = "$OUTPUT_ROOT/prank-codon-masked";
if(!-d $OUTPUT_MASKED_ALIGNMENTS){
	make_path($OUTPUT_MASKED_ALIGNMENTS) or die "Cannot create output directory";
}

### 
my $BASE_DIRECTORY = "$OUTPUT_GUIDANCE_PRANK_CODON";
opendir(BASEDIR, $BASE_DIRECTORY) or die $!;





###########################
####### MAIN SCRIPT #######
###########################
print STDERR "Masking PRANK multiple sequence alignments based on GUIDANCE results...\n\n";
my $count = 0;
my @removedSequenceSetsDueToLowGuidanceSequenceScore = ();

### loop over all guidance results directories in the base directory
while (my $basedir_file = readdir(BASEDIR)) {
	next unless(-d "$BASE_DIRECTORY/$basedir_file"); # only loop over directories
	next unless($basedir_file =~ /^(ENSG\d+)__cds$/); # only loop over relevant directories
	
	### PRINT COUNTER
	$count++;
	print STDERR $count . " MSAs done\r" if (($count % 100) == 0);


	### SET UP DIRECTORY
	my $current_ensembl_id = $1;
	# print $current_ensembl_id . "\n";

	my $CURRENT_GUIDANCE_DIRECTORY = "$BASE_DIRECTORY/$basedir_file";


	### CHECKS
	eval { checkPrankStdoutLog("$CURRENT_GUIDANCE_DIRECTORY" . "/" . "MSA.PRANK.aln.std") };
	if ($@){
		print $@;
		my $cmd = "mv $CURRENT_GUIDANCE_DIRECTORY $OUTPUT_GUIDANCE_PRANK_CODON_ERRORS/";
		print "\tMoving files to error directory: $cmd\n";
		`$cmd`;
		next;
	}


	### CHECKS
	eval { checkEndsOk("$CURRENT_GUIDANCE_DIRECTORY" . "/" . "ENDS_OK") };
	if ($@){
		print $@;
		my $cmd = "mv $CURRENT_GUIDANCE_DIRECTORY $OUTPUT_GUIDANCE_PRANK_CODON_ERRORS/";
		print "\tMoving files to error directory: $cmd\n";
		`$cmd`;
		next;
	}


	### GET ALIGNMENT
	my $alignment_base = undef;
	eval { $alignment_base = readBaseAlignment("$CURRENT_GUIDANCE_DIRECTORY" . "/" . "MSA.PRANK.aln.With_Names") };
	die if ($@);

	# my $AlignIO_out = Bio::AlignIO->new(-fh     => \*STDOUT ,
	# 									-format => 'fasta');
	# $AlignIO_out->write_aln($alignment_base);


	### GET ID MAPPING TABLES
	
	##
	#get mapping table between original sequence IDs and number IDs by GUIDANCE
	##
	my %sequenceIDMapping = (); # key: numbered ID by GUIDANCE, value: original sequence ID
	eval { %sequenceIDMapping = getSequenceIDMapping("$CURRENT_GUIDANCE_DIRECTORY" . "/" . "Seqs.Codes") };
	die if ($@);
	# printHash(\%sequenceIDMapping);

	##
	#get mapping table between alignment sequence IDs and the row number/order in the alignment (requried for residue-score mapping, see below)
	##
	my %alignmentIDRowMapping = getAlignmentIDRowMapping($alignment_base); # key: row number in alignment, value: original sequence ID at that row
	# printHash(\%alignmentIDRowMapping);


	### GET MASKING INFORMATION
	
	##
	#A) sequence-based
	##
	my $sequence_score_threshold = 0.6;
	my @sequencesToRemove = undef;

	eval { @sequencesToRemove = getSequencesToRemove("$CURRENT_GUIDANCE_DIRECTORY" . "/" . "MSA.PRANK.Guidance_res_pair_seq.scr_with_Names", $sequence_score_threshold) };
	die if ($@);
	# printArray(@sequencesToRemove);
	
	##
	#B) column-based
	##
	my $column_score_threshold = 0.93;
	my @columnsToRemove = undef;

	eval { @columnsToRemove = getColumnsToRemove("$CURRENT_GUIDANCE_DIRECTORY" . "/" . "MSA.PRANK.Guidance_res_pair_col.scr", $column_score_threshold) };
	die if ($@);
	# printArray(@columnsToRemove);

	##
	#C) residue-based
	##
	
	#NOTE: '-nan' values correspond to either a gap, or to a residue that needs to be removed, so before masking the residue, check if it's actually a gap, in which case keep it as a gap, or if it's a residue, in which case mask it
	my $residue_score_threshold = 0.9;
	my %residuesToRemove = (); # key: sequence ID, value: array of residues to be removed

	eval { %residuesToRemove = getResiduesToRemove("$CURRENT_GUIDANCE_DIRECTORY" . "/" . "MSA.PRANK.Guidance_res_pair_res.scr", $residue_score_threshold, \%alignmentIDRowMapping) };
	die if ($@);
	# printHashOfArrays(\%residuesToRemove);


	### MASK ALIGNMENTS
	my $GUIDANCE_MASK_CHARACTER = "n"; # reserved DNA chars: ATCG N
	
	##
	#A) sequence-based
	##
	# I want to keep all nine sequences; if any should be removed according to GUIDANCE, then throw away the whole set of orthologous sequences
	if(scalar @sequencesToRemove != 0){
		print STDERR "Removing sequence set $current_ensembl_id: one or more sequences have GUIDANCE score lower than $sequence_score_threshold\n";
		push(@removedSequenceSetsDueToLowGuidanceSequenceScore, $current_ensembl_id);
		next;
	}

	##
	#B) column-based
	##
	my $alignment_masked_columns = maskAlignmentColumns($alignment_base, \@columnsToRemove, $GUIDANCE_MASK_CHARACTER);
	
	##
	#C) residue-based
	##

	#NOTE: '-nan' values correspond to either a gap, or to a residue that needs to be removed, so before masking the residue, check if it's actually a gap, in which case keep it as a gap, or if it's a residue, in which case mask it
	my $alignment_masked_columns_and_residues = maskAlignmentResidues($alignment_masked_columns, \%residuesToRemove, $GUIDANCE_MASK_CHARACTER);



	### WRITE ALIGNMENTS
	writeAlignment($OUTPUT_MASKED_ALIGNMENTS . "/" . $current_ensembl_id . "__cds.prank-codon.aln.fa", $alignment_base);
	# writeAlignment($OUTPUT_MASKED_ALIGNMENTS . "/" . $current_ensembl_id . "__cds.prank-codon-guidance-masked-columns.aln.fa", $alignment_masked_columns);
	writeAlignment($OUTPUT_MASKED_ALIGNMENTS . "/" . $current_ensembl_id . "__cds.prank-codon-guidance-masked.aln.fa", $alignment_masked_columns_and_residues);


	### COPY OTHER FILES
	# my $cmd = ();
	# $cmd = "cp $CURRENT_GUIDANCE_DIRECTORY/MSA.PRANK.Guidance_res_pair_res.html $OUTPUT_MASKED_ALIGNMENTS/$current_ensembl_id" . "__cds.MSA.PRANK.Guidance_res_pair_res.html";
	# `$cmd`;

	# $cmd = "cp $OUTPUT_ROOT/$OUTPUT_CDS/$current_ensembl_id" . "__cds.prank-codon.aln.fa.best.fas $OUTPUT_MASKED_ALIGNMENTS/"; # the PRANK alignment in cds directory, not via GUIDANCE (don't use this one for masking)
	# `$cmd`;


	# closedir CURGUIDANCEDIR;
}
print STDERR $count . " MSAs done\n";
closedir BASEDIR;


### FINISH
print STDERR "\nRemoved " . scalar(@removedSequenceSetsDueToLowGuidanceSequenceScore) . " sequence sets because one or more sequences have GUIDANCE scores lower than the threshold\n";

my $removedSequenceSetsDueToLowGuidanceSequenceScoreFile = "removedSequenceSetsDueToLowGuidanceSequenceScore.txt";
open(REMOVEDSEQFILE, ">$removedSequenceSetsDueToLowGuidanceSequenceScoreFile");
print STDERR "\tWriting sequence sets to file $removedSequenceSetsDueToLowGuidanceSequenceScoreFile\n";
foreach(@removedSequenceSetsDueToLowGuidanceSequenceScore){
	print REMOVEDSEQFILE $_ . "\n";
}
close REMOVEDSEQFILE;

print STDERR "\nFINISHED\n";






#########################
####### FUNCTIONS #######
#########################
###########################################################################
sub writeAlignment{
	my ($alignment_file, $alignment) = @_;

	my $AlignIO_out = Bio::AlignIO->new(-file   => ">$alignment_file", #-fh     => \*STDOUT ,
	 									-format => 'fasta');
	$AlignIO_out->write_aln($alignment);
	# $AlignIO_out->write_aln($alignment->slice(1,1000,1));
}

###########################################################################
sub maskAlignmentResidues{
	my $alignment = shift @_;
	my %residuesToRemove = %{ shift @_ };
	my $mask_character = shift @_;

	# based on http://doc.bioperl.org/bioperl-live/Bio/SimpleAlign.html#CODE89
	my $alignment_masked_residues = Bio::SimpleAlign->new();
	$alignment_masked_residues->id($alignment->id);

	foreach my $seq ( $alignment->each_seq() ) {
		my $new_seq = Bio::LocatableSeq->new(   -id => $seq->id,
												-alphabet => $seq->alphabet,
												-verbose => $alignment->verbose);
		my $masked_seq_string = $seq->seq;

		foreach my $seq_id_of_residuesToRemove (keys %residuesToRemove){
			if($seq_id_of_residuesToRemove eq $seq->id){
				my @residues = @{ $residuesToRemove{$seq_id_of_residuesToRemove} };
				# print "Masking " . scalar @residues . " residues from $seq_id_of_residuesToRemove...\n";
				# print "\t[ " . join(",", @residues) . " ]\n";
				# separator3();
				
				foreach my $res (@residues){
					my $start = $res;
					my $end = $res;

					# convert from 1-based alignment coords!
					my $masked_seq_substring = substr($masked_seq_string, $start - 1, $end - $start + 1);

					#NOTE: '-nan' values correspond to either a gap, or to a residue that needs to be removed, so before masking the residue, check if it's actually a gap (or an undetermined character N from sequencing), in which case keep it as a gap, or if it's a residue, in which case mask it
					$masked_seq_substring =~ s/[^-N]/$mask_character/g; # replace by mask char, except when it's a gap

					# stitch together new sequence string
					$masked_seq_string = substr($masked_seq_string,0,$start-1) . $masked_seq_substring . substr($masked_seq_string,$end);
				}				
			}
		}

		$new_seq->seq($masked_seq_string);
		$alignment_masked_residues->add_seq($new_seq);
	}
	
	$alignment_masked_residues->set_displayname_flat(); # Function  : Makes all the sequences be displayed as just their name, not name/start-end
	return $alignment_masked_residues;
}

###########################################################################
sub maskAlignmentColumns{
	my $alignment = shift @_;
	my @columnsToRemove = @{ shift @_ };
	my $mask_character = shift @_;

	foreach my $col (@columnsToRemove){
		# print "Masking column #$col...\n";

		# print $Bio::LocatableSeq::GAP_SYMBOLS . "\n";
		# print $Bio::LocatableSeq::FRAMESHIFT_SYMBOLS . "\n";
		my $alignment_masked_columns = $alignment->mask_columns($col,$col,$mask_character);
		$alignment_masked_columns->set_displayname_flat(); # Function  : Makes all the sequences be displayed as just their name, not name/start-end

		$alignment = $alignment_masked_columns;
	}
	
	# my $AlignIO_out = Bio::AlignIO->new(-fh     => \*STDOUT ,
	#  									-format => 'fasta');
	# $AlignIO_out->write_aln($alignment->slice(1,50,1));
	
	return $alignment;
}

###########################################################################
sub getResiduesToRemove{
	my $file = shift @_;
	my $threshold = shift @_;
	my %alignmentIDRowMapping = %{shift @_};
	my %residuesToRemove = (); # key: sequence ID, value: array of residues to be removed

	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#read GUIDANCE sequence scores
	open(FILE, "<$file");
	while(my $line = <FILE>){
		next if($line =~ /^#/);
		chomp $line;
		
		if($line =~ /^\s*(\d+)\s+(\d+)\s+(.+)$/){
			my $seqColNumber = $1;
			
			my $seqRowNumber = $2;
			my $seqIDOriginal = $alignmentIDRowMapping{$seqRowNumber};

			my $residueScore = $3;


			#add residue to residues that are to be removed if it scores below the threshold, or is -nan
			#NOTE: '-nan' values correspond to either a gap, or to a residue that needs to be removed, so before masking the residue, check if it's actually a gap, in which case keep it as a gap, or if it's a residue, in which case mask it
			if($residueScore eq '-nan'){
				# push(@residuesToRemove, $1);
				push( @{ $residuesToRemove{$seqIDOriginal} }, $seqColNumber);
			} elsif($residueScore < $threshold){
				push( @{ $residuesToRemove{$seqIDOriginal} }, $seqColNumber);
			}

		} else {
			die "Could not parse line #$. \'$line\' in file $file\n";
		}
		
	}
	close FILE;

    return %residuesToRemove;
}

###########################################################################
sub getColumnsToRemove{
	my $file = shift @_;
	my $threshold = shift @_;
	my @columnsToRemove = ();

	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#read GUIDANCE sequence scores
	open(FILE, "<$file");
	while(my $line = <FILE>){
		next if($line =~ /^#/);
		chomp $line;
		
		if($line =~ /^\s*(\d+)\s+(.+)$/){

			#add column to columns that are to be removed if it scores below the threshold
			if($2 < $threshold){
				push(@columnsToRemove, $1);
			}

		} else {
			die "Could not parse line #$. \'$line\' in file $file\n";
		}
		
	}
	close FILE;

    return @columnsToRemove;
}

###########################################################################
sub getSequencesToRemove{
	my $file = shift @_;
	my $threshold = shift @_;
	my @sequencesToRemove = ();

	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#read GUIDANCE sequence scores
	open(FILE, "<$file");
	while(my $line = <FILE>){
		next if($line =~ /^SEQUENCE_NAME/ or $line =~ /^#END/);
		chomp $line;
		my @F = split /\t/, $line;
		
		# add sequence to sequences that are to be removed if it scores below the threshold
		if($F[1] < $threshold){
			push(@sequencesToRemove, $F[0]);
		}
	}
	close FILE;

    return @sequencesToRemove;
}


###########################################################################
sub getAlignmentIDRowMapping{
	my $alignment = shift @_;
	my %alignmentIDRowMapping = (); # key: row number in alignment, value: original sequence ID at that row

	my $row_number = 1;
	foreach my $seq ( $alignment->each_seq() ) {
		$alignmentIDRowMapping{$row_number} = $seq->display_id;
		$row_number += 1;
	}

    return %alignmentIDRowMapping;
}


###########################################################################
sub getSequenceIDMapping{
	my $file = shift @_;
	my %sequenceIDMapping = (); # key: numbered ID by GUIDANCE, value: original sequence ID

	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#read GUIDANCE sequence scores
	open(FILE, "<$file");
	while(my $line = <FILE>){
		chomp $line;
		
		my @F = split /\t/, $line;
		$F[1] =~ s/_placeholder_to_prevent_PRANK_error//;

		$sequenceIDMapping{$F[1]} = $F[0];		
	}
	close FILE;

    return %sequenceIDMapping;
}

###########################################################################
sub readBaseAlignment{
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
sub checkEndsOk{
	my $file = shift @_;
	
	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#if file is not empty (0 bytes, http://perldoc.perl.org/functions/-X.html) then GUIDANCE did not finish without error
	if(! -z $file){
		die "GUIDANCE did not \'end OK\': see $file";
	}
}

###########################################################################
sub checkPrankStdoutLog{
	my $file = shift @_;
	
	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	# test PRANK file
	local $/=undef;
	open(FILE, "<$file");
	my $filestring = <FILE>;
	close FILE;

	#file has to contain this pattern, otherwise something went wrong with PRANK
	if( $filestring !~ /Writing\s+- alignment to \'.+?\/MSA\.PRANK\.aln\.best\.fas\'/ ){
		die "PRANK problem: see $file";
	}
}

###########################################################################
sub checkIfExists{
	my $file = shift @_;
	if(! -f $file){
		die "File $file does not exist";
	}
}

