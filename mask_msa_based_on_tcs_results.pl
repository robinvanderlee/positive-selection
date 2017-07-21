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

require('functions.pl');





### output directory for sequence files
my $OUTPUT_ROOT = "sequences";
my $OUTPUT_CDS = "cds";
my $OUTPUT_TCS_PRANK_CODON = "$OUTPUT_ROOT/tcs-prank-codon";
my $OUTPUT_MASKED_ALIGNMENTS = "$OUTPUT_ROOT/prank-codon-masked";


### 
my $BASE_DIRECTORY = "$OUTPUT_TCS_PRANK_CODON";
opendir(BASEDIR, $BASE_DIRECTORY) or die $!;





###########################
####### MAIN SCRIPT #######
###########################
print STDERR "Masking PRANK multiple sequence alignments based on TCS results...\n\n";
print STDERR "NOTE THAT THE TCS RESULTS ARE ON AMINO ACID LEVEL, AND ALIGNMENT MAPPING NEEDS TO BE DONE ON CODON LEVEL...!!\n\n";
my $count = 0;
my @removedSequenceSetsDueToLowTCSSequenceScore = ();

### loop over all tcs results file in the base directory
while (my $basedir_file = readdir(BASEDIR)) {
	next unless(-f "$BASE_DIRECTORY/$basedir_file"); # only loop over files
	
	# only loop over relevant files, i.e. the TCS ASCII output files
	if($basedir_file =~ /^(ENSG\d+)__cds\.prank\-codon\.aln\.translated\.score_ascii$/){

		### PRINT COUNTER
		$count++;
		print STDERR $count . " MSAs done\r" if (($count % 100) == 0);

		### SET UP DIRECTORY
		my $current_ensembl_id = $1;
		# print $current_ensembl_id . "\n";
		


		### GET ALIGNMENT
		my $alignment_base = undef;
		eval { $alignment_base = readBaseAlignment("$OUTPUT_MASKED_ALIGNMENTS" . "/" . $current_ensembl_id . "__cds.prank-codon-guidance-masked.aln.fa") }; # GUIDANCE MASKED ALIGNMENT
		# eval { $alignment_base = readBaseAlignment("$OUTPUT_MASKED_ALIGNMENTS" . "/" . $current_ensembl_id . "__cds.prank-codon.aln.fa") }; # NON-MASKED PRANK-CODON ALIGNMENT
		die if ($@);

		# my $AlignIO_out = Bio::AlignIO->new(-fh     => \*STDOUT ,
		# 									-format => 'fasta');
		# $AlignIO_out->write_aln($alignment_base);


		### READ CONFIDENCE SCORES
		my @returns = ();
		eval { @returns = readTCSConfidenceScoresAscii("$BASE_DIRECTORY" . "/" . $basedir_file) };
		die if ($@);

		my %sequencesTCSConfidenceScoresAscii = %{ $returns[0] }; # keys: sequence IDs, values: TCS scores
		my $columnsTCSConfidenceScoresAscii = $returns[1]; # string of column scores
		my %residuesTCSConfidenceScoresAscii = %{ $returns[2] }; # keys: sequence IDs, strings of residue scores

		# printHash(\%sequencesTCSConfidenceScoresAscii);
		# separator2();
		# print $columnsTCSConfidenceScoresAscii . "\n";
		# separator2();
		# printHash(\%residuesTCSConfidenceScoresAscii);
		# separator2();


		### GET MASKING INFORMATION
		##
		#A) sequence-based
		##
		my $sequence_score_threshold = 50;
		my @sequencesToRemove = getSequencesToRemove(\%sequencesTCSConfidenceScoresAscii, $sequence_score_threshold);
		# printArray(@sequencesToRemove);
		
		##
		#B) column-based
		##
		my $column_score_threshold = 4;
		my @AAcolumnsToRemove = getColumnsToRemove($columnsTCSConfidenceScoresAscii, $column_score_threshold);
		# printArray(@AAcolumnsToRemove);

		##
		#C) residue-based
		##
		my $residue_score_threshold = 4;
		my %AAresiduesToRemove = getResiduesToRemove(\%residuesTCSConfidenceScoresAscii, $residue_score_threshold); # key: sequence ID, value: array of residues to be removed
		# printHashOfArrays(\%AAresiduesToRemove);


		### MASK ALIGNMENTS
		my $TCS_MASK_CHARACTER = "n"; # reserved DNA chars: ATCG N
		
		##
		#A) sequence-based
		##
		# I want to keep all nine sequences; if any should be removed according to TCS, then throw away the whole set of orthologous sequences
		if(scalar @sequencesToRemove != 0){
			print STDERR "Removing sequence set $current_ensembl_id: one or more sequences have TCS score lower than $sequence_score_threshold\n";
			push(@removedSequenceSetsDueToLowTCSSequenceScore, $current_ensembl_id);
			next;
		}

		##
		#B) column-based
		##
		my $alignment_masked_columns = maskAlignmentColumns($alignment_base, \@AAcolumnsToRemove, $TCS_MASK_CHARACTER);
			
		##
		#C) residue-based
		##
		my $alignment_masked_columns_and_residues = maskAlignmentResidues($alignment_masked_columns, \%AAresiduesToRemove, $TCS_MASK_CHARACTER);


		### WRITE ALIGNMENTS
		# writeAlignment($OUTPUT_MASKED_ALIGNMENTS . "/" . $current_ensembl_id . "__cds.prank-codon.aln.fa", $alignment_base);
		# writeAlignment($OUTPUT_MASKED_ALIGNMENTS . "/" . $current_ensembl_id . "__cds.prank-codon-guidance-masked-tcs-masked-columns.aln.fa", $alignment_masked_columns);
		writeAlignment($OUTPUT_MASKED_ALIGNMENTS . "/" . $current_ensembl_id . "__cds.prank-codon-guidance-tcs-masked.aln.fa", $alignment_masked_columns_and_residues);


		### COPY OTHER FILES
		# my $cmd = ();
		# $cmd = "cp $OUTPUT_ROOT/$OUTPUT_CDS/$current_ensembl_id" . "__cds.prank-codon.aln.fa.best.fas $OUTPUT_MASKED_ALIGNMENTS/"; # the PRANK alignment in cds directory, not via GUIDANCE (don't use this one for masking)
		# `$cmd`;

	}
}
print STDERR $count . " MSAs done\n";
closedir BASEDIR;


### FINISH
print STDERR "\nRemoved " . scalar(@removedSequenceSetsDueToLowTCSSequenceScore) . " sequence sets because one or more sequences have TCS scores lower than the threshold\n";

my $removedSequenceSetsDueToLowTCSSequenceScoreFile = "removedSequenceSetsDueToLowTCSSequenceScore.txt";
open(REMOVEDSEQFILE, ">$removedSequenceSetsDueToLowTCSSequenceScoreFile");
print STDERR "\tWriting sequence sets to file $removedSequenceSetsDueToLowTCSSequenceScoreFile\n";
foreach(@removedSequenceSetsDueToLowTCSSequenceScore){
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
	my %AAresiduesToRemove = %{ shift @_ };
	my $mask_character = shift @_;

	# based on http://doc.bioperl.org/bioperl-live/Bio/SimpleAlign.html#CODE89
	my $alignment_masked_residues = Bio::SimpleAlign->new();
	$alignment_masked_residues->id($alignment->id);

	foreach my $seq ( $alignment->each_seq() ) {
		my $new_seq = Bio::LocatableSeq->new(   -id => $seq->id,
												-alphabet => $seq->alphabet,
												-verbose => $alignment->verbose);
		my $masked_seq_string = $seq->seq;

		foreach my $seq_id_of_residuesToRemove (keys %AAresiduesToRemove){
			if($seq_id_of_residuesToRemove eq $seq->id){
				my @AAresidues = @{ $AAresiduesToRemove{$seq_id_of_residuesToRemove} };
				# print "Masking " . scalar @AAresidues . " residues from $seq_id_of_residuesToRemove...\n";
				# print "\t[ " . join(",", @AAresidues) . " ]\n";
				# separator3();
				
				foreach my $AAcol (@AAresidues){
					#	AA
					# 	1	2	3	4	5
					# 	Q	E	R	S	T
					#
					#	CODON
					# 	123	456	789	012	345
					# 	AGT	GCA	TTC	CCC	AGG
					#
					#	((a - 1) * 3) + 1
					my $CODONcol = (($AAcol - 1) * 3) + 1;

					my $start = $CODONcol;
					my $end = $CODONcol + 2;

					# convert from 1-based alignment coords!
					my $masked_seq_substring = substr($masked_seq_string, $start - 1, $end - $start + 1);

					#before masking the residue, check if it's actually a gap (or an undetermined character N from sequencing), in which case keep it as a gap, or if it's a residue, in which case mask it
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
	my @AAcolumnsToRemove = @{ shift @_ };
	my $mask_character = shift @_;

	foreach my $AAcol (@AAcolumnsToRemove){
		# print "Masking column #$col...\n";

		#	AA
		# 	1	2	3	4	5
		# 	Q	E	R	S	T
		#
		#	CODON
		# 	123	456	789	012	345
		# 	AGT	GCA	TTC	CCC	AGG
		#
		#	((a - 1) * 3) + 1
		my $CODONcol = (($AAcol - 1) * 3) + 1;

		# print $Bio::LocatableSeq::GAP_SYMBOLS . "\n";
		# print $Bio::LocatableSeq::FRAMESHIFT_SYMBOLS . "\n";
		my $alignment_masked_columns = $alignment->mask_columns($CODONcol,$CODONcol+2,$mask_character);
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
	my %residuesTCSConfidenceScoresAscii = %{ shift @_ }; # keys: sequence IDs, strings of residue scores
	my $threshold = shift @_;
	my %residuesToRemove = (); # key: sequence ID, value: array of residues to be removed

	foreach my $seqID (keys %residuesTCSConfidenceScoresAscii){
		my $residueScoreString = $residuesTCSConfidenceScoresAscii{$seqID};

		my @F = split //, $residueScoreString;
		for(my $i = 0; $i < (scalar @F); $i++){

			# add residue to residues that are to be removed if it scores below the threshold
			if($F[$i] ne '-' && $F[$i] < $threshold){
				push( @{ $residuesToRemove{$seqID} }, $i + 1); # @F is 0-based
				# print $F[$i] . "\t" . ($i + 1) . "\n"; 
			}
		}
	}

    return %residuesToRemove;
}

###########################################################################
sub getColumnsToRemove{
	my $columnsTCSConfidenceScoresAscii = shift @_; # string of column scores
	my $threshold = shift @_;
	my @columnsToRemove = ();

	my @F = split //, $columnsTCSConfidenceScoresAscii;
	for(my $i = 0; $i < (scalar @F); $i++){
		#add column to columns that are to be removed if it scores below the threshold
		if($F[$i] < $threshold){
			push(@columnsToRemove, $i + 1); # @F is 0-based
			# print $F[$i] . "\t" . ($i + 1) . "\n"; 
		}
	}

    return @columnsToRemove;
}

###########################################################################
sub getSequencesToRemove{
	my %sequencesTCSConfidenceScoresAscii = %{ shift @_ }; # keys: sequence IDs, values: TCS scores
	my $threshold = shift @_;
	my @sequencesToRemove = ();

	foreach(keys %sequencesTCSConfidenceScoresAscii){
		# add sequence to sequences that are to be removed if it scores below the threshold
		if($sequencesTCSConfidenceScoresAscii{$_} < $threshold){
			push(@sequencesToRemove, $_);
		}
	}

    return @sequencesToRemove;
}

###########################################################################
sub readTCSConfidenceScoresAscii{
	my $file = shift @_;

	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);


	### READ TCS CONFIDENCE SCORES
	open(FILE, "<$file");
	local $/ = "\n\n"; # separate file into blocks
	chomp (my @blocks = <FILE>);
	close FILE;



	### SEQUENCE INFO BLOCK
	#
	# T-COFFEE, Version_11.00.61eb9e4 (2015-04-14 13:25:34 - Revision 61eb9e4 - Build 495)
	# Cedric Notredame 
	# CPU TIME:0 sec.
	# SCORE=977
	# *
	#  BAD AVG GOOD
	# *
	# 9606__Homo_sapiens__ENSG00000197993__ENST00000355265__ENSP00000347409                       :  98
	# 9598__Pan_troglodytes__ENSPTRG00000024026__ENSPTRT00000049635__ENSPTRP00000043404           :  97
	# 9595__Gorilla_gorilla_gorilla__ENSGGOG00000007998__ENSGGOT00000008041__ENSGGOP00000007829   :  98
	# 9601__Pongo_abelii__ENSPPYG00000018107__ENSPPYT00000021109__ENSPPYP00000020313              :  98
	# 61853__Nomascus_leucogenys__ENSNLEG00000002907__ENSNLET00000003824__ENSNLEP00000003645      :  98
	# 9544__Macaca_mulatta__ENSMMUG00000006571__ENSMMUT00000009202__ENSMMUP00000008647            :  96
	# 9555__Papio_anubis__ENSPANG00000021000__ENSPANT00000014598__ENSPANP00000018713              :  95
	# 60711__Chlorocebus_sabaeus__ENSCSAG00000007359__ENSCSAT00000005370__ENSCSAP00000003592      :  98
	# 9483__Callithrix_jacchus__ENSCJAG00000009471__ENSCJAT00000018445__ENSCJAP00000017440        :  98
	# cons                                                                                        :  97
	#
	my $sequence_score_info = shift @blocks;
	my @sequence_score_info_lines = split /\n/, $sequence_score_info;

	my %sequencesTCSConfidenceScoresAscii = (); # keys: sequence IDs, values: TCS scores
	foreach my $line (@sequence_score_info_lines){
		next unless $line =~ /^(\S+)\s+\:\s+(\d+)$/;
		next if $line =~ /^cons/;

		$sequencesTCSConfidenceScoresAscii{$1} = $2;
	}


	### RESIDUE AND COLUMN INFO BLOCKS
	#
	# 9606__Homo_sapiens__ENSG00000197993__ENST00000355265__ENSP00000347409                       8999999999989999999999999998899988-------9999999
	# 9598__Pan_troglodytes__ENSPTRG00000024026__ENSPTRT00000049635__ENSPTRP00000043404           899999999998999999999999999-------11111129999999
	# 9595__Gorilla_gorilla_gorilla__ENSGGOG00000007998__ENSGGOT00000008041__ENSGGOP00000007829   899999999998999999999999999888888811111228899999
	# 9601__Pongo_abelii__ENSPPYG00000018107__ENSPPYT00000021109__ENSPPYP00000020313              89999999999-7899999999999998899988-------9999999
	# 61853__Nomascus_leucogenys__ENSNLEG00000002907__ENSNLET00000003824__ENSNLEP00000003645      --------------99999999999998899988-------9999999
	# 9544__Macaca_mulatta__ENSMMUG00000006571__ENSMMUT00000009202__ENSMMUP00000008647            -789999999989999999999999998899888-------9999999
	# 9555__Papio_anubis__ENSPANG00000021000__ENSPANT00000014598__ENSPANP00000018713              8999999999999999999999999998899988-------9999999
	# 60711__Chlorocebus_sabaeus__ENSCSAG00000007359__ENSCSAT00000005370__ENSCSAP00000003592      8999999999989999999999999998899988-------9999999
	# 9483__Callithrix_jacchus__ENSCJAG00000009471__ENSCJAT00000018445__ENSCJAP00000017440        ------------------------------------------------
	# cons                                                                                        889999999998999999999999999889988811111119999999
	#
	my $columnsTCSConfidenceScoresAscii = "";
	my %residuesTCSConfidenceScoresAscii = ();

	foreach my $block (@blocks){
		next if($block =~ /^\s*$/);
		my @residues_and_column_score_info_lines = split /\n/, $block;
		
		foreach my $line (@residues_and_column_score_info_lines){
			next if($line =~ /^\s*$/);
			
			## COLUMNS
			if($line =~ /^cons\s+(\d+)/){
				$columnsTCSConfidenceScoresAscii .= $1;

			## RESIDUES
			} elsif($line =~ /^(\S+)\s+([\d\-]+)$/){
				$residuesTCSConfidenceScoresAscii{$1} .= $2;

			##
			} else {
				die "Could not parse file $file, problem with block:\n$block\n\n";
			}
		}
	}

    return (\%sequencesTCSConfidenceScoresAscii, $columnsTCSConfidenceScoresAscii, \%residuesTCSConfidenceScoresAscii);
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
sub checkIfExists{
	my $file = shift @_;
	if(! -f $file){
		die "File $file does not exist";
	}
}

