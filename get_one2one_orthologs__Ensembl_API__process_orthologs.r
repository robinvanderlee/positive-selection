#####################################
## Robin van der Lee               ##
## robinvanderlee AT gmail DOT com ##
############################################################################################################
## Genome-scale detection of positive selection in 9 primates predicts human-virus evolutionary conflicts ##
## Robin van der Lee, Laurens Wiel, Teunis J.P. van Dam, Martijn A. Huynen                                ##
############################################################################################################

### INPUT FILE
#head -10 one2one_orthologs_primates_20150303-1449.txt 
# ENSG00000228913	
# ENSG00000275008	
# ENSG00000180383	pan_troglodytes	ENSPTRG00000013353	ENSPTRP00000022893	nomascus_leucogenys	ENSNLEG00000009576	ENSNLEP00000011652	chlorocebus_sabaeus	ENSCSAG00000012789	ENSCSAP00000008967	papio_anubis	ENSPANG00000016042	ENSPANP00000016887	macaca_mulatta	ENSMMUG00000005116	ENSMMUP00000006779	callithrix_jacchus	ENSCJAG00000020809	ENSCJAP00000038711	
# ENSG00000277657	
# ENSG00000162444	pongo_abelii	ENSPPYG00000001904	ENSPPYP00000002207	macaca_mulatta	ENSMMUG00000002219	ENSMMUP00000002963	chlorocebus_sabaeus	ENSCSAG00000000775	ENSCSAP00000017231	papio_anubis	ENSPANG00000015220	ENSPANP00000004180	callithrix_jacchus	ENSCJAG00000005470	ENSCJAP00000010001	nomascus_leucogenys	ENSNLEG00000009377	ENSNLEP00000011417	pan_troglodytes	ENSPTRG00000000125	ENSPTRP00000000237	gorilla_gorilla	ENSGGOG00000028335	ENSGGOP00000018605	
# ENSG00000165583	gorilla_gorilla	ENSGGOG00000010688	ENSGGOP00000010422	pongo_abelii	ENSPPYG00000020365	ENSPPYP00000022796	
# ENSG00000104671	pan_troglodytes	ENSPTRG00000020134	ENSPTRP00000044768	gorilla_gorilla	ENSGGOG00000013074	ENSGGOP00000012754	pongo_abelii	ENSPPYG00000018486	ENSPPYP00000020722	nomascus_leucogenys	ENSNLEG00000011906	ENSNLEP00000014459	macaca_mulatta	ENSMMUG00000013735	ENSMMUP00000018058	papio_anubis	ENSPANG00000008817	ENSPANP00000002958	chlorocebus_sabaeus	ENSCSAG00000015839	ENSCSAP00000011920	callithrix_jacchus	ENSCJAG00000012338	ENSCJAP00000022619	
# ENSG00000101400	pan_troglodytes	ENSPTRG00000013400	ENSPTRP00000022970	gorilla_gorilla	ENSGGOG00000008508	ENSGGOP00000008318	pongo_abelii	ENSPPYG00000010924	ENSPPYP00000012200	nomascus_leucogenys	ENSNLEG00000009970	ENSNLEP00000012154	papio_anubis	ENSPANG00000016612	ENSPANP00000016827	chlorocebus_sabaeus	ENSCSAG00000012284	ENSCSAP00000008479	macaca_mulatta	ENSMMUG00000012719	ENSMMUP00000016694	callithrix_jacchus	ENSCJAG00000018516	ENSCJAP00000034358	
# ENSG00000228691	
# ENSG00000196368	pan_troglodytes	ENSPTRG00000021910	ENSPTRP00000037622	gorilla_gorilla	ENSGGOG00000015735	ENSGGOP00000015349	nomascus_leucogenys	ENSNLEG00000008174	ENSNLEP00000009962	macaca_mulatta	ENSMMUG00000032455	ENSMMUP00000012438	papio_anubis	ENSPANG00000006680	ENSPANP00000009391
# ENSG00000228913    
# ENSG00000275008	
# ENSG00000180383	pan_troglodytes	ENSPTRG00000013353	ENSPTRP00000022893	nomascus_leucogenys	ENSNLEG00000009576	ENSNLEP00000011652	chlorocebus_sabaeus	ENSCSAG00000012789	ENSCSAP00000008967	papio_anubis	ENSPANG00000016042	ENSPANP00000016887	macaca_mulattENSMMUG00000005116	ENSMMUP00000006779	callithrix_jacchus	ENSCJAG00000020809	ENSCJAP00000038711	
orthologs <- read.table(unz("ensembl78_one2one_orthologs_primates_20150303-1449.txt.zip", "one2one_orthologs_primates_20150303-1449.txt"), header = F, sep = "\t", quote = "", na.strings = "", fill = T)
dim(orthologs)
head(orthologs)
summary(orthologs)
table(table(orthologs[,1]))

### remove last column, which is always empty
orthologs <- orthologs[,-26]

### discard rows where one or more species do not have an ortholog listed
# we search for orthologs to human genes (column 1) in 8 species
# if a species is found, there are 3 columns: species_name, ensembl gene id, ensembl protein id
num.species.required <- 8
num.cols.required <- num.species.required * 3 + 1


keep.discard <- apply(orthologs, MARGIN = 1, function(x){
    if(length(which(!is.na(x))) == num.cols.required){
        return("keep");
    } else {
        return("discard");
    }
})
keep.discard
table(keep.discard)

orthologs.all.species <- orthologs[keep.discard == "keep",]

dim(orthologs.all.species)
table(table(orthologs.all.species[,1])) # 11170 human genes with 1:1:1:1.. orthologs in all 9 primates

apply(orthologs.all.species, 2, function(x) table(table(x)))


### sort on human ID
orthologs.all.species <- orthologs.all.species[order(orthologs.all.species[,1]),]
head(orthologs.all.species)


### order columns
# biomart approach had this order:
# Ensembl.Gene.ID    Chimpanzee.Ensembl.Gene.ID	Gibbon.Ensembl.Gene.ID	Gorilla.Ensembl.Gene.ID	Marmoset.Ensembl.Gene.ID	Olive.baboon.Ensembl.Gene.ID	Orangutan.Ensembl.Gene.ID	Vervet.AGM.Ensembl.Gene.ID	Macaque.Ensembl.Gene.ID
table(as.vector(as.matrix(orthologs.all.species[,seq(2,25,3)])))
write.table(t(t(names(table(as.vector(as.matrix(orthologs.all.species[,seq(2,25,3)])))))), quote = F, row.names = F, col.names = F)

species.order <- c("pan_troglodytes",
                   "nomascus_leucogenys",
                   "gorilla_gorilla",
                   "callithrix_jacchus",
                   "papio_anubis",
                   "pongo_abelii",
                   "chlorocebus_sabaeus",
                   "macaca_mulatta")
species.order


get.ordered.orthologs <- function(x) {
    x <- as.vector(as.matrix(x))
    
    # exclude human element
    human.id <- x[1]
    x <- x[-1]
    
    split.factor <- sort(rep(1:(length(x)/3),3)) # groups for splitting
    split.list <- split(x, split.factor)
    order.hash <- unlist(lapply(split.list, function(y){  return(which(species.order == y[1]))  })) # check which position each list element is supposed to have
    sorted.hash <- sort(order.hash) # sorted hash, where the names now contain the correct order of split.list elements
    split.list.sorted <- split.list[as.integer(names(sorted.hash))] # sorted list
    
    c(human.id, as.vector(unlist(lapply(split.list.sorted, function(z){ return(z[2]) } ))))
}

orthologs.all.species.ordered.frame <- matrix(nrow = nrow(orthologs.all.species), ncol = 9)
for(i in 1:nrow(orthologs.all.species)){
    row <- orthologs.all.species[i,]
    orthologs.all.species.ordered.frame[i,] <- get.ordered.orthologs(row)
}

orthologs.all.species.ordered.frame <- as.data.frame(orthologs.all.species.ordered.frame)
names(orthologs.all.species.ordered.frame) <- c("homo_sapiens", species.order)
dim(orthologs.all.species.ordered.frame)
head(orthologs.all.species.ordered.frame)


### write tables
write.table(orthologs.all.species.ordered.frame, file = "ensembl_API_one2one_ortholog_clusters.txt", quote = F, sep = "\t", row.names = F, col.names = T)


