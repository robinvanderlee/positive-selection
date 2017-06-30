##############################
## Robin van der Lee        ##
## robinvanderlee@gmail.com ##
##############################

### INPUT FILES
#head -10 /Users/robinvanderlee/Research/analysis_positive_selection_primates/get_one2one_orthologs__biomart/mart_export__?.txt 
# ==> /Users/robinvanderlee/Research/analysis_positive_selection_primates/get_one2one_orthologs__biomart/mart_export__1.txt <==
# Ensembl Gene ID	Associated Gene Name	Chimpanzee Ensembl Gene ID	Canonical Protein or Transcript ID	Homology Type	Gibbon Ensembl Gene ID	Canonical Protein or Transcript ID	Homology Type	Gorilla Ensembl Gene ID	Canonical Protein or Transcript ID	Homology Type	Marmoset Ensembl Gene ID	Canonical Protein or Transcript ID	Homology Type	Olive baboon Ensembl Gene ID	Canonical Protein or Transcript ID	Homology Type	Orangutan Ensembl Gene ID	Canonical Protein or Transcript ID	Homology Type
# ENSG00000228913	UBD														
# ENSG00000275008	KIR2DL3														
# ENSG00000180383	DEFB124	ENSPTRG00000013353	ENSP00000326309	ortholog_one2one	ENSNLEG00000009576	ENSP00000326309ortholog_one2one				ENSCJAG00000020809	ENSP00000326309	ortholog_one2one	ENSPANG00000016042	ENSP00000326309	ortholog_one2one			
# ENSG00000277657	KIR2DS3														
# ENSG00000162444	RBP7	ENSPTRG00000000125	ENSP00000294435	ortholog_one2one	ENSNLEG00000009377	ENSP00000294435ortholog_one2one	ENSGGOG00000028335	ENSP00000294435	ortholog_one2one	ENSCJAG00000005470	ENSP00000294435	ortholog_one2one	ENSPANG00000015220	ENSP00000294435	ortholog_one2one	ENSPPYG00000001904	ENSP00000294435	ortholog_one2one
# ENSG00000165583	SSX5							ENSGGOG00000010688	ENSP00000312415	ortholog_one2onENSCJAG00000014776	ENSP00000312415	ortholog_one2many	ENSPANG00000006212	ENSP00000312415	ortholog_one2many	ENSPPYG00000020365	ENSP00000312415	ortholog_one2one
# ENSG00000104671	DCTN6	ENSPTRG00000020134	ENSP00000221114	ortholog_one2one	ENSNLEG00000011906	ENSP00000221114ortholog_one2one	ENSGGOG00000013074	ENSP00000221114	ortholog_one2one	ENSCJAG00000012338	ENSP00000221114	ortholog_one2one	ENSPANG00000008817	ENSP00000221114	ortholog_one2one	ENSPPYG00000018486	ENSP00000221114	ortholog_one2one
# ENSG00000101400	SNTA1	ENSPTRG00000013400	ENSP00000217381	ortholog_one2one	ENSNLEG00000009970	ENSP00000217381ortholog_one2one	ENSGGOG00000008508	ENSP00000217381	ortholog_one2one	ENSCJAG00000018516	ENSP00000217381	ortholog_one2one	ENSPANG00000016612	ENSP00000217381	ortholog_one2one	ENSPPYG00000010924	ENSP00000217381	ortholog_one2one
# ENSG00000228691	NEU1														

# ==> /Users/robinvanderlee/Research/analysis_positive_selection_primates/get_one2one_orthologs__biomart/mart_export__2.txt <==
# Ensembl Gene ID	Associated Gene Name	Vervet-AGM Ensembl Gene ID	Canonical Protein or Transcript ID	Homology Type	Macaque Ensembl Gene ID	Canonical Protein or Transcript ID	Homology Type
# ENSG00000228913	UBD						
# ENSG00000275008	KIR2DL3						
# ENSG00000180383	DEFB124	ENSCSAG00000012789	ENSP00000326309	ortholog_one2one	ENSMMUG00000005116	ENSP00000326309ortholog_one2one
# ENSG00000277657	KIR2DS3						
# ENSG00000162444	RBP7	ENSCSAG00000000775	ENSP00000294435	ortholog_one2one	ENSMMUG00000002219	ENSP00000294435ortholog_one2one
# ENSG00000165583	SSX5	ENSCSAG00000013354	ENSP00000312415	ortholog_one2many	ENSMMUG00000008046	ENSP00000312415ortholog_one2many
# ENSG00000104671	DCTN6	ENSCSAG00000015839	ENSP00000221114	ortholog_one2one	ENSMMUG00000013735	ENSP00000221114ortholog_one2one
# ENSG00000101400	SNTA1	ENSCSAG00000012284	ENSP00000217381	ortholog_one2one	ENSMMUG00000012719	ENSP00000217381ortholog_one2one
# ENSG00000228691	NEU1



### read data
biomart1 <- read.delim(unz("mart_export__1.txt.zip", "mart_export__1.txt"))
biomart2 <- read.delim(unz("mart_export__2.txt.zip", "mart_export__2.txt"))
dim(biomart1)
head(biomart1)
dim(biomart2)
head(biomart2)

### merge lists
biomart <- merge(biomart1,biomart2,all = T,by = c("Ensembl.Gene.ID", "Associated.Gene.Name"))
t(t(names(biomart)))
names(biomart)[grep(pattern = "Ensembl.Gene.ID", names(biomart))]
head(biomart,1)
dim(biomart)

table(table(biomart$Ensembl.Gene.ID))

### check distribution of orthology types
names(biomart)[grep(pattern = "Homology.Type", names(biomart))]
head(biomart[,grep(pattern = "Homology.Type", names(biomart))])
table(as.vector(as.matrix(biomart[,grep(pattern = "Homology.Type", names(biomart))])))

### discard rows that contain many2many and one2many ortholog relationships, so that we are left with 1:1:1:1..etc orthologous groups only
keep.discard <- apply(biomart, MARGIN = 1, function(x){
                                                if("ortholog_many2many" %in% x[grep(pattern = "Homology.Type", names(biomart))]){
                                                    return("discard");
                                                } else if("ortholog_one2many" %in% x[grep(pattern = "Homology.Type", names(biomart))]){
                                                    return("discard");
                                                } else {
                                                    return("keep");
                                                }
                                            })
biomart.one2one <- biomart[keep.discard == "keep",]

dim(biomart.one2one)
head(biomart.one2one)
table(table(biomart.one2one$Ensembl.Gene.ID)) # 3798 human genes have many2many or one2many relationships

table(!(levels(biomart.one2one$Ensembl.Gene.ID) %in% biomart.one2one$Ensembl.Gene.ID))
length(levels(biomart.one2one$Ensembl.Gene.ID)[!(levels(biomart.one2one$Ensembl.Gene.ID) %in% biomart.one2one$Ensembl.Gene.ID)])
biomart.NOTone2one <- levels(biomart.one2one$Ensembl.Gene.ID)[!(levels(biomart.one2one$Ensembl.Gene.ID) %in% biomart.one2one$Ensembl.Gene.ID)]

### discard rows where one or more species do not have an ortholog listed
keep.discard <- apply(biomart.one2one, MARGIN = 1, function(x){
    if("" %in% x[grep(pattern = "Homology.Type", names(biomart))]){
        return("discard");
    } else {
        return("keep");
    }
})
table(keep.discard) # 7015 human genes for which in one or multiple other primates no ortholog was annotated
biomart.one2one.all.species <- biomart.one2one[keep.discard == "keep",] 

dim(biomart.one2one.all.species)
table(table(biomart.one2one.all.species$Ensembl.Gene.ID)) # 11170 human genes with 1:1:1:1.. orthologs in all 9 primates

apply(biomart.one2one.all.species, 2, function(x) table(table(x)))

### remove columns Canonical.Protein.or.Transcript.ID
names(biomart.one2one.all.species)[grep(pattern = "Canonical.Protein.or.Transcript.ID", names(biomart.one2one.all.species))]
head(biomart.one2one.all.species[,names(biomart)[grep(pattern = "Canonical.Protein.or.Transcript.ID", names(biomart.one2one.all.species))]])

same.different <- apply(biomart.one2one.all.species, MARGIN = 1, function(x){ # check if all of these columns are indeed the same
    if(length(unique(x[grep(pattern = "Canonical.Protein.or.Transcript.ID", names(biomart.one2one.all.species))])) == 1 ){
        return("same");
    } else {
        return("different");
    }
})
table(same.different)

### write tables
write.table(biomart.one2one.all.species[,names(biomart.one2one.all.species)[grep(pattern = "Ensembl.Gene.ID", names(biomart.one2one.all.species))]], file = "biomart_one2one_ortholog_clusters.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(biomart.one2one.all.species[,c(2,grep(pattern = "Ensembl.Gene.ID", names(biomart.one2one.all.species)))], file = "biomart_one2one_ortholog_clusters_with_gene_names.txt", quote = F, sep = "\t", row.names = F, col.names = T)

