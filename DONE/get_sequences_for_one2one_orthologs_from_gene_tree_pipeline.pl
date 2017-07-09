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
use File::Path qw( make_path );
use Data::Dumper;
use Bio::SeqIO;
use Bio::AlignIO;

require('functions.pl');
require('Ensembl_API_subroutines__RvdL/get_genetree.pl');


### output directory for sequence files
my $OUTPUT_ROOT = "sequences";
my $OUTPUT_PROTEIN = "protein";
my $OUTPUT_CDS = "cds";

my @dirs = ("$OUTPUT_ROOT/$OUTPUT_PROTEIN", "$OUTPUT_ROOT/$OUTPUT_CDS");
foreach my $dir(@dirs){
  if(!-d $dir){
    make_path($dir) or die "Cannot create output directories";
  }
}





###################################################################################################
### read one2one orthologs
# Ensembl.Gene.ID Chimpanzee.Ensembl.Gene.ID  Gibbon.Ensembl.Gene.ID  Gorilla.Ensembl.Gene.ID Marmoset.Ensembl.Gene.ID  Olive.baboon.Ensembl.Gene.ID  Orangutan.Ensembl.Gene.ID Vervet.AGM.Ensembl.Gene.ID  Macaque.Ensembl.Gene.ID
# ENSG00000000003 ENSPTRG00000022087  ENSNLEG00000007505  ENSGGOG00000010003  ENSCJAG00000004288  ENSPANG00000018367  ENSPPYG00000020538  ENSCSAG00000010385  ENSMMUG00000008209
# ENSG00000000005 ENSPTRG00000022086  ENSNLEG00000007504  ENSGGOG00000017007  ENSCJAG00000004270  ENSPANG00000016641  ENSPPYG00000020537  ENSCSAG00000010392  ENSMMUG00000008212
# ENSG00000000419 ENSPTRG00000013625  ENSNLEG00000006998  ENSGGOG00000008993  ENSCJAG00000000551  ENSPANG00000000450  ENSPPYG00000011131  ENSCSAG00000015216  ENSMMUG00000002759
# ...
my $ortholog_file = "biomart_one2one_ortholog_clusters.txt";
open(IN, "<$ortholog_file") or die "Cannot open $ortholog_file\n";
while(<IN>){};
my $total_genes = $. - 1; # minus header
close IN;

print "Need to get sequences for a total of " . $total_genes . " orthologs.\n";
# print "Already obtained sequences for " . scalar @genes_done_previously . " orthologs in previous runs.\n";
separator();






###################################################################################################
### which genes does this script need to do?

my @genes_to_do = ();
if(defined $ARGV[0] and $ARGV[0] eq "-p"){ # passes as arugment?
  ### parallel

  if(defined $ARGV[1] and $ARGV[1] eq "-f"){
    my $parallel_genes_file = $ARGV[2];
    open(PAR, "<$parallel_genes_file") or die "Parallel requested, but no file with genes provided!\n";
    while(<PAR>){
      chomp;
      push(@genes_to_do,$_);
    }
    close PAR;

  } else {
    die "Parallel requested, but no file with genes provided!\n";
  }  

} else {
  # ### test genes
  # # MAVS
  # push(@genes_to_do, "ENSG00000088888");

  push(@genes_to_do, "ENSG00000019549");

  # print join(", ", @genes_to_do) . "\n";
  # separator();
}
my %genes_to_do_hash = map { $_ => 1 } @genes_to_do; # make unique

my $total_gene_to_do = scalar @genes_to_do;
print "Getting sequences of orthologs for " . $total_gene_to_do . " genes in this run...\n";
separator();





### get registry
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;

print "Connecting to Ensembl API registry...\n";
printf( "The API version used is %s\n", software_version() );
separator();
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(
  -host=>'ensembldb.ensembl.org',
  -user=>'anonymous', 
  );
$reg->set_reconnect_when_lost(); # should prevent errors like: DBD::mysql::st execute failed: Lost connection to MySQL server during query at $PATH_TO_Ensembl_API/ensembl/modules/Bio/EnsEMBL/DBSQL/TranscriptAdaptor.pm line 238, <IN> line 25.
# TODO CHECK AND PRINT ENSEMBL VERSION



### get adaptors
print "Obtaining Ensembl API adaptors...\n";
separator();
my $genemember_adaptor = $reg->get_adaptor('Multi', 'compara', 'GeneMember');
my $comparaDBA = $reg->get_DBAdaptor('compara', 'compara');
my $genetree_adaptor = $comparaDBA->get_GeneTreeAdaptor;
my $gene_adaptor = $reg->get_adaptor( "human", "core", "gene" );
print "genemember_adaptor schema version: " . $genemember_adaptor->schema_version() . "\n";
print "genetree_adaptor schema version: " . $genetree_adaptor->schema_version() . "\n";
print "gene_adaptor schema version: " . $gene_adaptor->schema_version() . "\n";



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


### which genes were already done?
# my $genes_done_file = "$OUTPUT_ROOT/genes_done.txt";
# my @genes_done_previously = ();
# if (-f $genes_done_file) {
#   open(DONEPREV, "<$genes_done_file");
#   while(<DONEPREV>){
#     chomp;
#     push(@genes_done_previously, $_);
#   }
#   close DONEPREV;
# }
# my %genes_done_previously_hash = map { $_ => 1 } @genes_done_previously;










###################################################################################################

### MAIN procedure
### for each human gene (1 per row), fetch the genetree objects, and from those:
###   1. fetch the sequences for the orthologs on the same line as the human gene 
###   2. fetch the compara alignment sequences for these orthologs
open(IN, "<$ortholog_file") or die "Cannot open $ortholog_file\n";
my @genes_done = ();
while (my $line = <IN>){
  next if $. == 1; # skip header
  chomp $line;

  my @one2one_ortholog_cluster = split /\t/, $line;
  my %one2one_ortholog_cluster_hash = map { $_ => 1 } @one2one_ortholog_cluster;
  
  my $human_gene_id = $one2one_ortholog_cluster[0];

  # next if exists $genes_done_previously_hash{$human_gene_id}; # already obtained the sequences for this gene in a previous run
  next unless exists $genes_to_do_hash{$human_gene_id};


  print "Getting GeneTree for $human_gene_id...\n";
  separator2();
  
  my $human_gene_id_genemember = $genemember_adaptor->fetch_by_stable_id($human_gene_id);
  die "Cannot fetch GeneMember for $human_gene_id\n" unless (defined $human_gene_id_genemember);

  my $genetree = get_genetree($genetree_adaptor, $human_gene_id_genemember);
  die "Cannot fetch GeneTree for $human_gene_id\n" unless (defined $human_gene_id_genemember);
  

  print "Getting GeneTree members...\n";
  separator2();
  my @genetree_members = @{ $genetree->get_all_Members };


  my @seq_types = ('protein','cds');
  my @sequence_objects_protein = ();
  my %sequence_objects_protein_ids_hash = ();
  my @sequence_objects_cds = ();

  print "Looping over GeneTree members...\n";
  separator2();

  ### loop over genetree_members to find those members that are in the ortholog list, and get their sequences
  while ( my $genetree_member = shift @genetree_members ){
    my $genetree_member_gene_member = $genetree_member->gene_member();
    my $genetree_member_transcript_member = $genetree_member->get_Transcript();

    my $taxon_object = $genetree_member_gene_member->taxon();
    my $current_taxon_id = $taxon_object->taxon_id;
    my $current_binomial = $taxon_object->binomial;
    $current_binomial =~ s/ /_/g;
    
    my $common_name = "";
    if(defined $taxon_object->common_name){
      $common_name = $taxon_object->common_name;
    }

    ### select only genetree_members matching the one2one ortholog selection
    if(exists $one2one_ortholog_cluster_hash{$genetree_member_gene_member->stable_id}){
      print $genetree_member->stable_id . " - ";
      print $genetree_member_transcript_member->stable_id . " - ";
      print $genetree_member_gene_member->stable_id . " - ";
      print $current_binomial . " - ";
      print $common_name . " - ";
      print $current_taxon_id;
      nl();


      # http://www.ensembl.org/info/docs/Doxygen/compara-api/classBio_1_1EnsEMBL_1_1Compara_1_1SeqMember.html#abcdd0ba1c3b5c7462a6c86a93b033dab
      #  public String Bio::EnsEMBL::Compara::SeqMember::other_sequence  (   ) 
      #  Arg [1]    : string $seq_type
      #  Arg [2]    : (opt) string $sequence
      #  Example    :
      # my $filtered_seq = $member->other_sequence('cds');
      #  Description: Get/Set the alternative sequence of type $seq_type for this member.
      #                Currently, proteins have 'cds', 'exon_bounded' and 'exon_cased'
      #                sequences, while RNAs have 'seq_with_flanking' sequences.
      #                The undef $seq_type maps to the default sequence ($member->sequence())
      #                If $sequence is set: store it in the database.
      #                'exon_cased' maps to the sequence string of this member with alternating upper
      #                and lower case corresponding to the translateable exons.
      #                'exon_bounded' maps to the sequence string of this member with exon boundaries
      #                denoted as O, B, or J depending on the phase (O=0, B=1, J=2)
      
      ### get individual protein and cDNA sequences
      foreach my $seq_type (@seq_types){
        my $alphabet = "";
        my $sequence = "";

        if($seq_type eq 'protein'){
          $alphabet = 'protein';
          
          # get the sequence
          $sequence = $genetree_member->sequence();
        }elsif($seq_type eq 'cds'){
          $alphabet = 'dna';

          # get the sequence
          $sequence = $genetree_member->other_sequence('cds');
        }

        # construct sequence ID
        my $seqname = $current_taxon_id . "__" .
                         $current_binomial . "__" .
                         $genetree_member_gene_member->stable_id . "__" .
                         $genetree_member_transcript_member->stable_id . "__" .
                         $genetree_member->stable_id;


        # construct bioperl sequence object
        my $sequence_object = Bio::Seq->new(
                                            -seq        => $sequence,
                                            -display_id => $seqname,
                                            -accession_number  => $current_taxon_id, # store taxon in accession_number
                                        );

        if($seq_type eq 'protein'){
          push(@sequence_objects_protein, $sequence_object);
        }elsif($seq_type eq 'cds'){
          push(@sequence_objects_cds, $sequence_object);
        }
      }

      $sequence_objects_protein_ids_hash{$genetree_member->stable_id} = 1; # build hash of ensembl protein IDs that are in my sequence object array

      separator2();
    }

  }

  ### write protein and cds sequence arrays to file
  # first order them on species
  my @sequence_objects_protein_ordered = order_sequences_by_taxon(\@sequence_objects_protein, \@species_of_interest_display_order);
  my @sequence_objects_cds_ordered = order_sequences_by_taxon(\@sequence_objects_cds, \@species_of_interest_display_order);

  foreach my $seq_type (@seq_types){
    my $filename = $human_gene_id . "__" . $seq_type . ".fa";
    my $seqout = Bio::SeqIO->new
                                (-file => ">$OUTPUT_ROOT/$seq_type/$filename",
                                 -format => 'fasta');

    if($seq_type eq 'protein'){
      $seqout->write_seq(@sequence_objects_protein_ordered);
    }elsif($seq_type eq 'cds'){
      $seqout->write_seq(@sequence_objects_cds_ordered);
    }
    print "File $filename has been generated\n";
  }





  ###### get Compara multiple alignments of protein and cDNA sequences in the genetree,
  ######     filter for species of interest, remove all-gap columns, sort sequences
  separator2();
  print "Getting Compara alignments...\n";
  separator2();

  foreach my $seq_type (@seq_types){
    my $seq_type_tmp = $seq_type;
    if($seq_type eq 'protein'){
      $seq_type_tmp = undef;
    }

    my $alignment = $genetree->get_SimpleAlign
                                       (-SEQ_TYPE => $seq_type_tmp,
                                       -ID_TYPE => 'STABLE',
                                       -STOP2X => 'false',
                                       -APPEND_TAXON_ID => 'true',
                                       -APPEND_SP_SHORT_NAME => 'true',
                                       -KEEP_GAPS => 'false');

    ### select specific sequences
    foreach my $alignment_sequence ($alignment->each_seq) {
      my @F = split /\_/, $alignment_sequence->display_id;
      my $ensemble_protein_id = $F[0];
      my $alignment_sequence_taxon = $F[1];

      # if(!exists $species_of_interest{$alignment_sequence_taxon}){
      #   $alignment->remove_seq($alignment_sequence);
      #   next;
      # }

      # check if the protein ID exists in the genetree sequences fetched above
      if(!exists $sequence_objects_protein_ids_hash{$ensemble_protein_id}){
        $alignment->remove_seq($alignment_sequence);
      }
    }

    ### no gaps-only columns
    $alignment = $alignment->remove_columns(['all_gaps_columns']);

    ### order alignment sequences on species
    my @sequence_objects_compara_alignment = ();
    foreach my $alignment_sequence ($alignment->each_seq) {
      my @F = split /\_/, $alignment_sequence->display_id;
      my $alignment_sequence_taxon = $F[1];

      my $sequence_object = Bio::Seq->new(
                                    -seq        => $alignment_sequence->seq,
                                    -display_id => $alignment_sequence->display_id,
                                    -accession_number  => $alignment_sequence_taxon, # store taxon in accession_number
                                );
      push(@sequence_objects_compara_alignment, $sequence_object)
    }
    my @sequence_objects_compara_alignment_ordered = order_sequences_by_taxon(\@sequence_objects_compara_alignment, \@species_of_interest_display_order);

    ### write alignment to file
    my $filename = $human_gene_id . "__" . $seq_type . ".compara.aln.fa";
    my $seqout = Bio::SeqIO->new
                                (-file => ">$OUTPUT_ROOT/$seq_type/$filename",
                                 -format => 'fasta');

    $seqout->write_seq(@sequence_objects_compara_alignment_ordered);
    print "File $filename has been generated\n";
  }

  push(@genes_done, $human_gene_id);
  separator3();
  print scalar @genes_done . " of $total_gene_to_do genes done...\n";
  separator3();
}
close IN;


### update file with the list of genes that was done:
# update_genes_done_file($genes_done_file, \@genes_done, \@genes_done_previously);


### disconnect registry
$reg->clear(); 


### finish
separator();
print "FINISHED\n";









###################################################################################################
sub order_sequences_by_taxon{
  # first argument: array of Bio::Seq objects, with -accession_number set as taxon
  # second argument: array of taxonomy numbers
  my $seq_objects_ref = shift @_;
  my @seq_objects = @{$seq_objects_ref};

  my $species_order_ref = shift @_;
  my @species_order = @{$species_order_ref};

  my @seq_objects_ordered = ();
  foreach my $species_order_taxon (@species_order){
    foreach my $seq_object (@seq_objects){
      if($species_order_taxon == $seq_object->accession_number){ # taxon number is stored in accession_number!
        push(@seq_objects_ordered, $seq_object);
      }
    }
  }

  return @seq_objects_ordered;
}

sub update_genes_done_file{
  my $genes_done_file = shift @_;
  my $genes_done_ref = shift @_;
  my $genes_done_previously_ref = shift @_;

  my @genes_done = @{$genes_done_ref};
  my @genes_done_previously = @{$genes_done_previously_ref};

  open(DONENEW, ">$genes_done_file");
  my @genes_done_total = (@genes_done_previously, @genes_done);
  my %genes_done_total_hash = map { $_ => 1 } @genes_done_total; # get uniques
  foreach(keys %genes_done_total_hash){
    print DONENEW "$_\n";
  }
  close DONENEW;
}
