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
require('functions.pl');




#####################
####### SETUP #######
#####################
### command line arguments
my $OUTPUT_SPECIFIC_DIRECTORY = "";
if(scalar @ARGV == 0){
	die "Please provide a codeml output directory\n";
} else {
	$OUTPUT_SPECIFIC_DIRECTORY = shift @ARGV;
}

### base output directory for codeml results directories
my $OUTPUT_ROOT = ".";

###
my $OUTPUT_DIRECTORY_FULL_PATH = "$OUTPUT_ROOT/$OUTPUT_SPECIFIC_DIRECTORY";
opendir(BASEDIR, $OUTPUT_DIRECTORY_FULL_PATH) or die $!;

### parse ensembl ID
$OUTPUT_DIRECTORY_FULL_PATH =~ /(ENSG\d+)/;
my $current_ensembl_id = $1;

### FILES WITH PARSED DATA
my $RESULTS_PARSED_OUTPUT_ROOT = "$OUTPUT_ROOT/codeml_results_parsed";
`mkdir $RESULTS_PARSED_OUTPUT_ROOT` unless -d $RESULTS_PARSED_OUTPUT_ROOT; # unless it already exists

# 1) general results on the alignment: e.g. kappa, omega, lnL, P values, etc.
my $PARSED_DATA_FILE_1__ALIGNMENT_CODEML_RESULTS = $RESULTS_PARSED_OUTPUT_ROOT . "/" . "$current_ensembl_id" . ".alignment_codeml_results";
open(PARSED_DATA_FILE_1, ">$PARSED_DATA_FILE_1__ALIGNMENT_CODEML_RESULTS");

# 2) residues with alternative model omegas, probabilities, etc.
my $PARSED_DATA_FILE_2__RESIDUES_CODEML_RESULTS = $RESULTS_PARSED_OUTPUT_ROOT . "/" . "$current_ensembl_id" . ".residues_codeml_results";
open(PARSED_DATA_FILE_2, ">$PARSED_DATA_FILE_2__RESIDUES_CODEML_RESULTS");




# R-E-S-U-L-T
# print PARSED_DATA_FILE_1 "alignment_ensembl_id\t" . $current_ensembl_id . "\n";




###########################
####### MAIN SCRIPT #######
###########################
print STDERR "Processing codeml results from directory $OUTPUT_DIRECTORY_FULL_PATH\n\n";
# print $current_ensembl_id . "\n";

### loop over all files in the codeml results directory
# but first I need to read the ctl file, otherwise I can't use the mapping stuff in the other functions below)
my %mapping_table__originalpos_char;
my %mapping_table__originalpos_gaplesspos;

while (my $basedir_file = readdir(BASEDIR)) {
	my $basedir_file_full_path = "$OUTPUT_DIRECTORY_FULL_PATH/$basedir_file";
	next unless(-f $basedir_file_full_path); # only loop over files

	if( $basedir_file_full_path =~ /\.ctl$/ ){
		eval {
			my $alignment_file = parseCtlFile($basedir_file_full_path);
			my $alignment_object = readAlignmentFile($alignment_file);

			`cp $alignment_file $RESULTS_PARSED_OUTPUT_ROOT`;

			# # see Bio::AlignIO::phylip
			# my $AlignIO_out = Bio::AlignIO->new(-fh     => \*STDOUT ,
			#  									-format 			=> 'phylip',
			#  									-idlength	 		=> 30,
			#  									# -longid		 		=> 1,
			#  									-flag_SI 			=> 1,
			#  									-interleaved 		=> 1,
			#  									-wrap_sequential	=> 1,
			#  									-idlinebreak 		=> 1);
			# $AlignIO_out->write_aln($alignment_object);

			my @mappingTableRefs = calculatePositionMappingTables($current_ensembl_id, $alignment_object);
			%mapping_table__originalpos_char = %{$mappingTableRefs[0]};
			%mapping_table__originalpos_gaplesspos = %{$mappingTableRefs[1]};
		};
		if ($@){
			print $@;
			die; # next;
		}
	}
}
closedir BASEDIR;

# then I can do the rest
opendir(BASEDIR, $OUTPUT_DIRECTORY_FULL_PATH) or die $!;
while (my $basedir_file = readdir(BASEDIR)) {
	my $basedir_file_full_path = "$OUTPUT_DIRECTORY_FULL_PATH/$basedir_file";
	next unless(-f $basedir_file_full_path); # only loop over files
	# next unless($basedir_file =~ /^(ENSG\d+)__cds$/); # only loop over relevant files
	# print $basedir_file . "\n";


	### PARSE & CHECKS & RESULTS
	if( $basedir_file_full_path =~ /\.screen_output$/ ){
		`cp $basedir_file_full_path $RESULTS_PARSED_OUTPUT_ROOT`;

		eval { parseScreenOutput($basedir_file_full_path) };
		if ($@){
			print $@;
			die; # next;
		}
	}

	if( $basedir_file_full_path =~ /\.codeml_outfile$/ ){
		`cp $basedir_file_full_path $RESULTS_PARSED_OUTPUT_ROOT`;

		eval { parseCodemlOutfile($basedir_file_full_path, $current_ensembl_id, \%mapping_table__originalpos_char, \%mapping_table__originalpos_gaplesspos) };
		if ($@){
			print $@;
			die; # next;
		}
	}

	if( $basedir_file_full_path =~ /rub$/ ){
		eval { parseRub($basedir_file_full_path, $current_ensembl_id) };
		if ($@){
			print $@;
			die; # next;
		}
	}
}

### FINISH
closedir BASEDIR;
close(PARSED_DATA_FILE_1);
close(PARSED_DATA_FILE_2);

# print STDERR "\nFINISHED\n";
print STDERR "=" x 100 . "\n";




#########################
####### FUNCTIONS #######
#########################
###########################################################################
sub parseBEB_results{
	my $file = shift @_;
	my $current_ensembl_id = shift @_;

	my $BEB_resultsString = shift @_;
	my @BEB_resultsStringLines = split "\n+", $BEB_resultsString;

	my %mapping_table__originalpos_char = %{shift @_};
	my %mapping_table__originalpos_gaplesspos = %{shift @_};
		
	# Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)
	# Positively selected sites (*: P>95%; **: P>99%)
	# (amino acids refer to 1st sequence: 9606__Homo_sapiens)
	#
	#             Pr(w>1)     post mean +- SE for w
	#
	#    101 N      0.620         2.040 +- 1.376
	#    145 T      0.821         2.610 +- 1.276
	#    182 A      0.693         2.249 +- 1.368
	#    196 P      0.924         2.871 +- 1.102
	my @residueNumber = ();
	my @residueType = ();
	my @residueProbabilityOmegaBiggerThanOne = ();
	my @residueProbabilityOmegaBiggerThanOneSignificance = ();
	my @residueOmegaBiggerThanOne = ();

	# my $P_CUTOFF = 0.9;
	
	# get reported residues (P > 0.5 are reported)
	foreach my $line (@BEB_resultsStringLines){
		if($line =~ /amino acids refer to 1st sequence/){
			if($line !~ /9606__Homo_sapiens/){
				die "Problem with sequence ordering in alignment: see $file";
			}
		}

		if($line =~ /^\s+(\d+)\s+(\w)\s+(\S+?)(\**)\s+(\S+)/){
			# print $line . "\n";
			
			# select residues based on P cutoff
			# if($3 > $P_CUTOFF){
				push(@residueNumber, $1);
				push(@residueType, $2);
				push(@residueProbabilityOmegaBiggerThanOne, $3);
				push(@residueProbabilityOmegaBiggerThanOneSignificance, $4);
				push(@residueOmegaBiggerThanOne, $5);
			# }
		}
	}

	# printArray(@residueNumber);
	# printArray(@residueType);

	# renumbering / position mapping
	my @gaplessResidueNumbers;
	for(my $i = 0; $i < scalar(@residueType); $i++){
		my $codemlResidueNumber = $residueNumber[$i];
		my $codemlResidueType = $residueType[$i];

		my $gaplessResidueNumber;
		if($mapping_table__originalpos_char{$codemlResidueNumber} eq $codemlResidueType){
			
			$gaplessResidueNumber = $mapping_table__originalpos_gaplesspos{$codemlResidueNumber};

			# print $codemlResidueNumber . "\n";
			# print $mapping_table__originalpos_char{$codemlResidueNumber} . "\n";
			# print $codemlResidueType . "\n";
			# print $gaplessResidueNumber . "\n";
			# print "\n";
		}
		push(@gaplessResidueNumbers, $gaplessResidueNumber);
	}


	# R-E-S-U-L-T
	for(my $i = 0; $i < scalar(@residueType); $i++){

		print PARSED_DATA_FILE_2	$current_ensembl_id . "\t" . 
									$residueType[$i] . "\t" .
									$residueNumber[$i] . "\t" .
									$gaplessResidueNumbers[$i] . "\t" .
									$residueProbabilityOmegaBiggerThanOne[$i] . "\t" .
									$residueProbabilityOmegaBiggerThanOneSignificance[$i] . "\t" .
									$residueOmegaBiggerThanOne[$i] . "\n";
	}
}

###########################################################################
sub parseKappaOmegaResults{
	my $kappaOmegaResultsString = shift @_;
	my $file = shift @_;
	my $current_ensembl_id = shift @_;
	
	my @kw_resultsStringLines = split "\n+", $kappaOmegaResultsString;
	
	my $kappa;
	my $omega;

	# ##OLD
	# # Detailed output identifying parameters
	# #
	# # kappa (ts/tv) =  5.17988
	# #
	# # Parameters in M8 (beta&w>1):
	# #   p0 =   0.97489  p =   0.25380 q =   1.64252
	# #  (p1 =   0.02511) w =   2.56165
	# foreach my $line (@kw_resultsStringLines){
	# 	# print $line . "\n";

	# 	if($line =~ /^kappa\s+\(ts\/tv\)\s+\=\s+(\S+)\s*$/){
	# 		$kappa = $1;
	# 	}
	# 	if($line =~ /w\s+\=\s+(\S+)\s*$/){
	# 		$omega = $1;
	# 	}
	# }

	#### NEW
	# M7vM8 ====================================== #
	# Detailed output identifying parameters

	# kappa (ts/tv) =  5.73281

	# Parameters in M8 (beta&w>1):
	#   p0 =   0.99986  p =   0.64638 q =  13.98859
	#  (p1 =   0.00014) w = 985.33511

	# dN/dS (w) for site classes (K=11)

	# p:   0.09999  0.09999  0.09999  0.09999  0.09999  0.09999  0.09999  0.09999  0.09999  0.09999  0.00014
	# w:   0.00060  0.00335  0.00764  0.01344  0.02100  0.03080  0.04379  0.06190  0.09015  0.15116 985.33609
	#
	# M1avM2a ====================================== #
	# Detailed output identifying parameters

	# kappa (ts/tv) =  3.84092

	# dN/dS (w) for site classes (K=3)

	# p:   0.76143  0.00000  0.23857
	# w:   0.00000  1.00000  1.17980
	#
	foreach my $line (@kw_resultsStringLines){
		# print $line . "\n";

		if($line =~ /^kappa\s+\(ts\/tv\)\s+\=\s+(\S+)\s*$/){
			$kappa = $1;
		}
		if($line =~ /^w\:\s+/){
			my @s = split /\s+/, $line;
			$omega = $s[$#s]; # get last element of split
		}
	}

	# check omega values
	if($omega == -1 || $omega == 89){
		die "codeml problem, indicated by omega value $omega: see $file";
	}

	# R-E-S-U-L-T
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" . "alternative_model_kappa\t" . $kappa . "\n";
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" . "alternative_model_omega\t" . $omega . "\n";
}

###########################################################################
sub parseLikelihoodsAndDoChi2{
	my @models = @{shift @_};
	my @likelihoodStrings = @{shift @_};
	my $file = shift @_;
	my $current_ensembl_id = shift @_;


	my @numberParameters = ();
	my @likelihoods = ();
	foreach my $match (@likelihoodStrings){

		if($match =~ /np\:\s+(\d+)\)\:\s+(-?[0-9]\d*(\.\d+)?)/){
			# number of parameters
			push(@numberParameters, $1);
			push(@likelihoods, $2);
		} else {
			die "Problem with the chi2 likelihood ratio test: see $file";
		}
	}
	# printArray(@likelihoodStrings);
	# printArray(@numberParameters);
	# printArray(@likelihoods);

	# do chi2 test // likelihood ratio test
	my $degreesOfFreedom = $numberParameters[1] - $numberParameters[0];
	my $LRT_statistic = -2 * ($likelihoods[0] - $likelihoods[1]); # see notes
	my $chi2_output = `chi2 $degreesOfFreedom $LRT_statistic`; # actually runs the PAML chi2 program
	
	my $prob1;
	my $prob2;
	if($chi2_output =~ /prob\s+\=\s+(\S+)\s+\=\s+(\S+)\s*$/){
		$prob1 = $1;
		$prob2 = $2;
	} elsif($chi2_output =~ /invalid/){
		# this happens when the LRT statistic is below smaller than zero (i.e. alternative model lnL is SMALLER than simple model lnL.. which should actually not be possible since the alternative model has more parameters and hence should always be able to fit the data better)
		$prob1 = "NA";
		$prob2 = "NA";
	} else {
		die "Problem with the chi2 likelihood ratio test: see $file";
	}

	# # print table
	# print 	"\t" .
	# 		"likelihood\t" .
	# 		"LRT_statistic [2x(lnL1 - lnL0)]\t" .
	# 		"parameters\t" .
	# 		"df\t" .
	# 		"P value" . "\n";

	# print 	$models[0] . "\t" .
	# 		$likelihoods[0] . "\t" .
	# 		$numberParameters[0] . "\t" .
	# 		"\t" .
	# 		"\t" .
	# 		"\t" . "\n";

	# print 	$models[1] . "\t" .
	# 		$likelihoods[1] . "\t" .
	# 		$numberParameters[1] . "\t" .
	# 		$LRT_statistic . "\t" .
	# 		$degreesOfFreedom . "\t" .
	# 		"$prob1 ($prob2)" . "\n";

	# R-E-S-U-L-T
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" . "null_model\t" . $models[0] . "\n";
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "null_model_lnL\t" . $likelihoods[0] . "\n";
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "null_model_np\t" . $numberParameters[0] . "\n";

	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "alternative_model\t" . $models[1] . "\n";
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "alternative_model_lnL\t" . $likelihoods[1] . "\n";
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "alternative_model_np\t" . $numberParameters[1] . "\n";
	
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "LRT_statistic\t" . $LRT_statistic . "\n";
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "LRT_degrees_of_freedom\t" . $degreesOfFreedom . "\n";
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "LRT_P_value_full\t" . $prob1 . "\n";
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "LRT_P_value_scientific\t" . $prob2 . "\n";
}

###########################################################################
sub parseCodemlOutfile{
	my $file = shift @_;
	my $current_ensembl_id = shift @_;
	my $oref = shift @_;
	my $gref = shift @_;

	my %mapping_table__originalpos_char = %{ $oref };
	my %mapping_table__originalpos_gaplesspos = %{ $gref };


	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#slurp file
	local $/=undef;
	open(FILE, "<$file");
	my $filestring = <FILE>;
	close FILE;

	## test if file contains this stuff
	my $countTimeUsed = () = $filestring =~ /Time used\:/gi;
	if( $countTimeUsed != 2 ){ # check if this pattern occurs twice
		die "codeml problem: see $file";
	}

	my $convergenceString = "check convergence\.\.";
	my @matches = ( $filestring =~ /$convergenceString/g );
	
	# R-E-S-U-L-T
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "codeml_outfile_convergence_warnings\t" . scalar(@matches) . "\n";

	my @models = ( $filestring =~ /Model \d\: (.+)\n/g );
	# printArray(@models);
	if( scalar(@models) != 2 ){ # check if this pattern occurs twice
		die "codeml problem: see $file";
	}

	## maximum likelihood estimations for the models
	my @likelihoodLines = ( $filestring =~ /(lnL\(ntime\:.+)\n/g );
	eval{ parseLikelihoodsAndDoChi2(\@models, \@likelihoodLines, $file, $current_ensembl_id) };
	die if($@);
	
	## get kappa and omega estimations for alternative model (the second occurence)
	my @kappaOmegaResultsStrings = ( $filestring =~ /(Detailed output identifying parameters.*?)dN \& dS for each branch/sg );
	eval{ parseKappaOmegaResults($kappaOmegaResultsStrings[1], $file, $current_ensembl_id) }; # only get the values for the second occurence (alternative model)
	die if($@);

	## get part of outfile with BEB results (positively selected residues)
	my $BEB_resultsString = "";
	if( $filestring =~ /(Bayes Empirical Bayes \(BEB\) analysis.*)The grid/s ){
		$BEB_resultsString = $1;
		eval{ parseBEB_results($file, $current_ensembl_id, $BEB_resultsString, \%mapping_table__originalpos_char, \%mapping_table__originalpos_gaplesspos) };
		die if($@);
	} else {
		die "codeml problem, no BEB results reported: see $file";
	}
}

###########################################################################
sub parseRub{
	my $file = shift @_;
	my $current_ensembl_id = shift @_;

	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#slurp file
	local $/=undef;
	open(FILE, "<$file");
	my $filestring = <FILE>;
	close FILE;

	#test if file contains this stuff
	my $convergenceString = "check convergence!";
	my @matches = ( $filestring =~ /$convergenceString/g );
	
	# R-E-S-U-L-T
	print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "rub_convergence_warnings\t" . scalar(@matches) . "\n";
}

###########################################################################
sub parseScreenOutput{
	my $file = shift @_;
	
	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#slurp file
	local $/=undef;
	open(FILE, "<$file");
	my $filestring = <FILE>;
	close FILE;

	#test if file contains this stuff
	my $countTimeUsed = () = $filestring =~ /Time used\:/gi;
	if( $countTimeUsed != 2 ){ # check if this pattern occurs twice
		die "codeml problem: see $file";
	}

	if( $filestring =~ /This is a rooted tree, without clock\.  Check\./ ){
        die "Rooted tree problem: see $file";
	}

	if( $filestring !~ /Branch lengths in tree are fixed\./ ){
		die "codeml problem: see $file";
	}
	
	my @matches = ( $filestring =~ /Model \d\: (\S+)/g );
	# printArray(@matches);
	if( scalar(@matches) != 2 ){ # check if this pattern occurs twice
		die "codeml problem: see $file";
	}
	
	if( $filestring !~ /BEBing \(dim \= 4\)\.  This may take several minutes\./ ){
		die "codeml problem: see $file";
	}
}

###########################################################################
sub calculatePositionMappingTables{
	my $current_ensembl_id = shift @_;
	my $alignment_object = shift @_;

	my %mapping_table__originalpos_char;
	my %mapping_table__originalpos_gaplesspos;

	foreach my $seq ( $alignment_object->each_seq() ){
		# get human sequence
		next unless $seq->id =~ /^9606__Homo_sapiens/;
		
		# print $seq->alphabet . "\n";
		# print $seq->seq . "\n";

		my $counter_original = 0;
		my $counter_gapless = 0;

		my @characters = split //, $seq->seq;
		foreach my $char (@characters){
			$counter_original = $counter_original + 1;

			if($char !~ /\-/){
				$counter_gapless = $counter_gapless + 1;
			}
			
			# print $char . "\t" . $counter_original . "\t" . $counter_gapless . "\n";

			$mapping_table__originalpos_char{$counter_original} = $char;
			$mapping_table__originalpos_gaplesspos{$counter_original} = $counter_gapless;
		}

		print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "length_alignment\t" . $counter_original . "\n";
		print PARSED_DATA_FILE_1 $current_ensembl_id . "\t" .  "length_human_sequence\t" . $counter_gapless . "\n";
	}

	return (\%mapping_table__originalpos_char, \%mapping_table__originalpos_gaplesspos);
}

###########################################################################
sub readAlignmentFile{
	my $file = shift @_;
	# print $file;

	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	#read alignment file
	my $alignment;
	my $AlignIO_in  = Bio::AlignIO->new(-file   => $file ,
                                	    -format => 'fasta');
    $alignment = $AlignIO_in->next_aln(); # Returns : a Bio::Align::AlignI compliant object
	$alignment->set_displayname_flat(); # Function  : Makes all the sequences be displayed as just their name, not name/start-end

    #throw error if file contains another alignment
    my $alignment2 = $AlignIO_in->next_aln();
    if(defined $alignment2){
    	die "$file contains multiple alignments";
    }
    
    return $alignment;
}

###########################################################################
sub parseCtlFile{
	my $file = shift @_;
	
	#stop function if file doesn't exists
	eval{ checkIfExists($file) };
	die if($@);

	open(FILE, "<$file");
	my $alignment_file_full_path;
	while(<FILE>){
		if(/seqfile\s+\=\s+(\S+)\s*/){
			$alignment_file_full_path = $1; #../../../../sequences/prank-codon-masked/ENSG00000122965__cds.prank-codon-guidance-tcs-masked-species-sorted.aln.phy
		}
	}
	close FILE;

	my @s = split /\//, $alignment_file_full_path;
	my $alignment_file = join("/", @s[2 .. $#s]);
	# my $alignment_file = join("/", @s[3 .. $#s]);
	# print $alignment_file;

	my $translated_alignment_file = $alignment_file;
	$translated_alignment_file =~ s/\.phy/.translated.fa/;
	
	return $translated_alignment_file;
}

###########################################################################
sub checkIfExists{
	my $file = shift @_;
	if(! -f $file){
		die "File $file does not exist";
	}
}

