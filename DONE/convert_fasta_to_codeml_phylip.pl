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
my $fasta_alignment_file = "";
if(scalar @ARGV == 0){
	die "Please provide a fasta file\n";
} else {
	$fasta_alignment_file = shift @ARGV;
	die "Please provide a fasta file\n" if $fasta_alignment_file !~ /fa[s]?$/;
}

### read fasta file
my $alignment_object = undef;
eval { $alignment_object = readFastaFile($fasta_alignment_file) };
die if ($@);

### check compatability with codeml
eval { checkCodemlCompatibility($fasta_alignment_file, $alignment_object) };
die if ($@);

### convert undetermined and masked amino acids to codeml ambiguity character: ?
my $codeml_ambiguity_char = "?";
# print STDERR "Converting undetermined and masked amino acids to codeml ambiguity character: $codeml_ambiguity_char ...\n";
my $alignment_codeml_ambiguity = convertToCodemlAmbiguity($fasta_alignment_file, $alignment_object, $codeml_ambiguity_char);

### write in phylip format
my $phylip_alignment_file = $fasta_alignment_file;
$phylip_alignment_file =~ s/\.fa[s]?/\.phy/;
writePhylipFile($phylip_alignment_file, $alignment_codeml_ambiguity);

# print STDERR "Finished!\n";





###########################################################################
sub writePhylipFile{
	my ($outfile, $alignment_object) = @_;

	# see Bio::AlignIO::phylip
	my $AlignIO_out = Bio::AlignIO->new(-file   			=> ">$outfile", #-fh     => \*STDOUT ,
	 									-format 			=> 'phylip',
	 									-idlength	 		=> 30,
	 									# -longid		 		=> 1,
	 									-flag_SI 			=> 1,
	 									-interleaved 		=> 1,
	 									-wrap_sequential	=> 1,
	 									-idlinebreak 		=> 1);
	$AlignIO_out->write_aln($alignment_object);
}

###########################################################################
sub convertToCodemlAmbiguity{
	my $file = shift @_;
	my $alignment = shift @_;
	my $codeml_ambiguity_char = shift @_;
	
	# based on http://doc.bioperl.org/bioperl-live/Bio/SimpleAlign.html#CODE89
	my $alignment_codeml_ambiguity = Bio::SimpleAlign->new();
	$alignment_codeml_ambiguity->id($alignment->id);

	foreach my $seq ( $alignment->each_seq() ) {
		my $new_seq = Bio::LocatableSeq->new(   -id => $seq->id,
												-alphabet => $seq->alphabet,
												-verbose => $alignment->verbose);
		my $ambiguity_seq_string = $seq->seq;
		
		$ambiguity_seq_string =~ s/[Nn]/$codeml_ambiguity_char/g; # replace undetermined and masked amino acids to codeml ambiguity character: ?

		$new_seq->seq($ambiguity_seq_string);
		$alignment_codeml_ambiguity->add_seq($new_seq);
	}
	
	$alignment_codeml_ambiguity->set_displayname_flat(); # Function  : Makes all the sequences be displayed as just their name, not name/start-end

	return $alignment_codeml_ambiguity;
}

###########################################################################
sub checkCodemlCompatibility{
	my $file = shift @_;
	my $alignment = shift @_;
	
	foreach my $seq ( $alignment->each_seq() ){
	
		#throw error if sequence name contains special characters
		if($seq->display_id =~ /[\/\\"',:#\(\)\$= ]/){
			die "Alignment in $file contains sequence names with special characters:\n" . "\t" . $seq->display_id . "\n";
		}


		#throw error if sequence contains stop codon
		my $translated_seq = $seq->translate(-complete_codons => 1);
		if($translated_seq->seq =~ /\*/){
			die "Alignment in $file contains sequence with stop codon:\n" . "\t" . $seq->display_id . "\n" . "\t" . $translated_seq->seq . "\n";
		}
	}

	#throw error if sequence contains special characters
	# print join("\n", sort $alignment->symbol_chars(1)) . "\n";
	if(join("", sort $alignment->symbol_chars(1)) =~ /[^ATCGNn-]/){
		die "Alignment in $file contains sequences with special characters: " . join("", sort $alignment->symbol_chars(1)) . "\n";
	}
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

	# remove the ensembl identifiers in the sequence IDs to make all .phy files compatbility with the species names in the same phylogenetic tree
# 	 9 1728 I
# 9483__Callithrix_jacchus__ENSC   ATGGACGATG ACTATGAAGC GTACCACAGT CTGTTCTTGT CGCTGCTCGG ACTCTGCCCG 
# 9606__Homo_sapiens__ENSG000001   ATGGACGACG ACTATGAAGC GTACCACAGT CTGTTCTTGT CGCTGCTCGG ACTCTGCCCG 
# 9598__Pan_troglodytes__ENSPTRG   ATGGACGACG ACTATGAAGC GTACCACAGT CTGTTCTTGT CGCTGCTCGG ACTCTGCCCG 
# 9595__Gorilla_gorilla_gorilla_   ATGGACGACG ACTATGAAGC GTACCACAGT CTGTTCTTGT CGCTGCTCGG ACTCTGCCCG 
# 9601__Pongo_abelii__ENSPPYG000   ATGGACGACG ACTATGAAGC GTACCACAGT CTGTTCTTGT CGCTGCTCGG ACTCTGCCCG 
# 61853__Nomascus_leucogenys__EN   ATGGACCACG ACTATGAAGC GTACCACAGT CTGTTCTTGT CGCTGCTCGG ACTCTGCCCG 
# 9544__Macaca_mulatta__ENSMMUG0   ATGGACGACG ACTATGAAGC GTACCACAGT CTGTTCTTGT CGCTGCTGGG ACTCTGCCCA 
# 9555__Papio_anubis__ENSPANG000   ATGGACGACG ACTATGAAGC GTACCACAGT CTGTTCTTGT CGCTGCTGGG ACTCTGCCCA 
# 60711__Chlorocebus_sabaeus__EN   ATGGACGACG ACTATGAAGC GTACCACAGT CTGTTCTTGT CGCTGCTGGG ACTCTGCCCA 
	foreach my $seq ( $alignment->each_seq() ){
		my $new_id = $1 if( $seq->display_id =~ /^(\d+__\S+?)__/ );
		$seq->display_id($new_id);
	}

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
