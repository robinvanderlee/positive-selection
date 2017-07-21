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
require('functions.pl');
require('Ensembl_API_subroutines__RvdL/get_human_protein_coding_genes.pl');
local $| = 1; # to allow print output to the command line to be updated


### get registry
use Bio::EnsEMBL::Registry;
my $reg = 'Bio::EnsEMBL::Registry';
$reg->load_registry_from_db(
  -host=>'ensembldb.ensembl.org',
  -user=>'anonymous',
);


### get adaptors
my $genomedb_adaptor = $reg->get_adaptor("Multi", "compara", "GenomeDB");
my $homology_adaptor = $reg->get_adaptor("Multi", "compara", "Homology");


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


### get all human protein coding genes
my $all = 1;
my @gene_ids = ();
if($all == 1){
  @gene_ids = get_human_protein_coding_genes($reg);
} else {
  # RIG-I, test gene
  @gene_ids = ("ENSG00000107201");
  # MAVS, test gene
  @gene_ids = ("ENSG00000088888", "ENSG00000107201");
  print join(", ", @gene_ids) . "\n";
}
separator();
my $total_genes = scalar @gene_ids;
print "Getting one2one orthologs for " . $total_genes . " genes...\n";
separator();


### for each species, fetch the GenomeDB oject
my %gdbs;
foreach my $species_taxon_num (keys %species_of_interest) {
  my $species_binomial = $species_of_interest{$species_taxon_num};
  $gdbs{$species_taxon_num} = $genomedb_adaptor->fetch_by_name_assembly($species_binomial);
}

# print name of the GenomeDB objects
print "Selected species:\n";
foreach my $this_genome_taxon_id (keys %gdbs){
  my $this_genome_db = $gdbs{$this_genome_taxon_id};
  print "\t" . $this_genome_db->name, "\n";
}
separator();





###################################################################################################
### Fetch homology object for a gene
# Homology Objects
# A Homology object represents either an orthologous or paralogous relationships between two members.
my $genemember_adaptor = $reg->get_adaptor('Multi', 'compara', 'GeneMember');


my $outfile = "one2one_orthologs_primates_" . getTimestampMin() . ".txt";
print "Writing one2one orthologs to file:\t$outfile\n\n";
open(OUT, ">$outfile") or die "Cannot open $outfile\n";


my $num_genes_done = 0;
select OUT;
### for each gene in the list, fetch all one2one orthologs in species of interest
foreach my $human_gene_id (@gene_ids){
  my $member = $genemember_adaptor->fetch_by_stable_id($human_gene_id);
  die "Cannot fetch GeneMember for $human_gene_id\n" unless (defined $member);

  print $member->stable_id . "\t";

  # get the homologies where the member is involved
  # That will return a reference to an array with all homologies (orthologues in
  # other species and paralogues in the same one)
  my $homologies = $homology_adaptor->fetch_all_by_Member($member);


  # Then for each homology, you can get all the Members implicated
  while ( my $homology = shift @{$homologies} ) {
    # http://www.ensembl.org/info/docs/Doxygen/compara-api/classBio_1_1EnsEMBL_1_1Compara_1_1Homology.html

    next unless $homology->description =~ /one2one/; # only select one2one orthologs!
    # print $homology->description . "\n";
    
    # Each homology relation has exactly 2 members, you should find there the initial member used as a query. The get_all_Members method returns an array of SeqMember objects. The SeqMember is actually an AlignedMember (for the underlying protein) and contains information about how this SeqMember has been aligned.
    my @homology_members = @{$homology->get_all_Members()};
    my $has_original_gene = 0;    
    foreach my $homology_member (@homology_members) {
      
      # each AlignedMember contains both the information on the SeqMember and in relation to the homology
      my $homology_member_gene_member = $homology_member->gene_member();

      my $id = $homology_member -> stable_id;
      my $geneid = $homology_member_gene_member -> stable_id;

      my $homology_member_genome_db = $homology_member->genome_db;
      die "Could not get GenomeDB for homology member $id\n" unless defined $homology_member_genome_db;
      my $homology_member_genome_db_taxon_id = $homology_member_genome_db->taxon_id;
      my $homology_member_genome_db_name = $homology_member_genome_db->name;

      if($geneid eq $human_gene_id){
        $has_original_gene = 1;
        die "Original gene should be human, but is " . $homology_member_genome_db_name . "\n" if $homology_member_genome_db_taxon_id != 9606;
      } else {
        # filter for species of interest
        foreach my $species_taxon_num (keys %species_of_interest){
          if($homology_member_genome_db_taxon_id eq $species_taxon_num){
            print $homology_member_genome_db_name . "\t";
            print $geneid . "\t";
            print $id . "\t";
            last;
          }
        }
      }

    }

    die "Original gene $human_gene_id was not part of the homology pair\n" if($has_original_gene == 0);
  }
  print "\n";

  $num_genes_done = $num_genes_done + 1;
  print STDOUT "$num_genes_done of $total_genes genes done...\r";
}
close OUT;
select STDOUT;
print "\n\n" . $outfile . " written\n";

separator();
print "FINISHED\n";

