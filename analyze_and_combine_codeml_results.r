#####################################
## Robin van der Lee               ##
## robinvanderlee AT gmail DOT com ##
############################################################################################################
## Genome-scale detection of positive selection in 9 primates predicts human-virus evolutionary conflicts ##
## Robin van der Lee, Laurens Wiel, Teunis J.P. van Dam, Martijn A. Huynen                                ##
############################################################################################################

#####################
####### SETUP #######
#####################
library("plyr")


#########################
####### FUNCTIONS #######
#########################
####### ALIGNMENT DATA #######
analyze.alignment_codeml_results <- function(codeml.parameter.combination){
    
    print("==========================================================================")    
        
    ### read data
    file.alignment_codeml_results <- paste("codeml_results_combined",
                                           "/",
                                           "codeml_results_base_files",
                                           "/",
                                           codeml.parameter.combination,
                                           ".alignment_codeml_results",
                                           sep = "")
    print(file.alignment_codeml_results)
    aln.data <- read.table(file = file.alignment_codeml_results, header = F, sep = "\t", colClasses = c("factor", "factor", "character"))
    
    names(aln.data) <- c("ensembl.id", "info.category", "info.value")
    head(aln.data)
    dim(aln.data)
    str(aln.data)
    
    
    ### split according to the type of information
    as.data.frame(levels(aln.data$info.category))
    #         levels(aln.data$info.category)
    # 1                     LRT_P_value_full
    # 2               LRT_P_value_scientific
    # 3               LRT_degrees_of_freedom
    # 4                        LRT_statistic
    # 5                    alternative_model
    # 6              alternative_model_kappa
    # 7                alternative_model_lnL
    # 8                 alternative_model_np
    # 9              alternative_model_omega
    # 10 codeml_outfile_convergence_warnings
    # 11                    length_alignment
    # 12               length_human_sequence
    # 13                          null_model
    # 14                      null_model_lnL
    # 15                       null_model_np
    # 16            rub_convergence_warnings
    
    lapply(split(aln.data, f = aln.data$info.category), nrow)
    lapply(split(aln.data, f = aln.data$info.category), head, 2)
    
    aln.data.split <- split(aln.data, f = aln.data$info.category)
    lapply(aln.data.split, function(x) class(x$info.value))
    
    
    ### analyze the new data frames as the correct data types
    
    # doesn't work # lapply(lapply(aln.data.split, function(x) do.call(data.frame, aln.data.split[[1]])), str)
    
    # doesn't work
    # sapply(names(aln.data.split), function(x){
    #         if(x %in% info.category.integer){
    #             tmp <- as.integer(aln.data.split[[x]]$info.value)
    #             aln.data.split[[x]]$info.value <- tmp
    #         } else {
    #             print("N")
    #         }
    #     })
    
    info.category.integer <- c("null_model_np", "alternative_model_np", "LRT_degrees_of_freedom", "codeml_outfile_convergence_warnings", "rub_convergence_warnings", "length_alignment", "length_human_sequence")
    info.category.character <- c("null_model", "alternative_model")
    info.category.numeric <- c("null_model_lnL", "alternative_model_lnL", "alternative_model_kappa", "alternative_model_omega", "LRT_statistic", "LRT_P_value_full", "LRT_P_value_scientific")
    
    
    ### some checks
    sapply(info.category.integer, function(x) table(as.integer(aln.data.split[[x]]$info.value)))
    sapply(info.category.character, function(x) table(as.character(aln.data.split[[x]]$info.value)))
    sapply(info.category.numeric, function(x) summary(as.numeric(aln.data.split[[x]]$info.value)))
    
    plot((aln.data.split[["LRT_P_value_full"]]$info.value), (aln.data.split[["LRT_P_value_scientific"]]$info.value))
    table((aln.data.split[["codeml_outfile_convergence_warnings"]]$info.value), (aln.data.split[["rub_convergence_warnings"]]$info.value))
    
    
    ### construct data frame with LRT info
    aln.data.split.LRT.list <-
        list(a = data.frame(ensembl.id = aln.data.split[["length_alignment"]][,"ensembl.id"],
                    length_alignment = as.numeric(aln.data.split[["length_alignment"]][,"info.value"])),
            b = data.frame(ensembl.id = aln.data.split[["length_human_sequence"]][,"ensembl.id"],
                    length_human_sequence = as.numeric(aln.data.split[["length_human_sequence"]][,"info.value"])),
            c = data.frame(ensembl.id = aln.data.split[["alternative_model_kappa"]][,"ensembl.id"],
                    alternative_model_kappa = as.numeric(aln.data.split[["alternative_model_kappa"]][,"info.value"])),
            d = data.frame(ensembl.id = aln.data.split[["alternative_model_omega"]][,"ensembl.id"],
                    alternative_model_omega = as.numeric(aln.data.split[["alternative_model_omega"]][,"info.value"])),
            e = data.frame(ensembl.id = aln.data.split[["null_model_lnL"]][,"ensembl.id"],
                    null_model_lnL = as.numeric(aln.data.split[["null_model_lnL"]][,"info.value"])),
            f = data.frame(ensembl.id = aln.data.split[["alternative_model_lnL"]][,"ensembl.id"],
                    alternative_model_lnL = as.numeric(aln.data.split[["alternative_model_lnL"]][,"info.value"])),
            g = data.frame(ensembl.id = aln.data.split[["LRT_statistic"]][,"ensembl.id"],
                    LRT_statistic = as.numeric(aln.data.split[["LRT_statistic"]][,"info.value"])),
            h = data.frame(ensembl.id = aln.data.split[["LRT_P_value_full"]][,"ensembl.id"],
                    LRT_P_value_full = as.numeric(aln.data.split[["LRT_P_value_full"]][,"info.value"])))
    lapply(aln.data.split.LRT.list, nrow)
    lapply(aln.data.split.LRT.list, head, 2)
    
    # combine all info into one data frame
    aln.data.LRT <- join_all(aln.data.split.LRT.list, by = "ensembl.id")
    head(aln.data.LRT)
    
    
    ### P value stuff: multiple testing correction etc.
    # for negative LRT test statistics, codeml would assign P value = NA, but since NA P values would not be taken into account by p.adjust, I converted NA P values to 1
    aln.data.LRT[which(is.na(aln.data.LRT[,"LRT_P_value_full"])),"LRT_P_value_full"] <- 1
    aln.data.LRT[,"LRT_P_value_log10"] <- log10(aln.data.LRT$LRT_P_value_full)
    aln.data.LRT[,"LRT_P_value_bonf"] <- p.adjust(aln.data.LRT$LRT_P_value_full, method = "bonf")
    aln.data.LRT[,"LRT_P_value_bh"] <- p.adjust(aln.data.LRT$LRT_P_value_full, method = "BH")
    
    head(aln.data.LRT)
    dim(aln.data.LRT)
    str(aln.data.LRT)
    
    
    ### distributions and some basic analysis
    summary(aln.data.LRT$alternative_model_kappa)
    hist(aln.data.LRT$alternative_model_kappa, breaks = 10000)
    hist(aln.data.LRT$alternative_model_kappa, xlim = c(0,10), breaks = 10000)
    summary(aln.data.LRT$alternative_model_omega)
    hist(aln.data.LRT$alternative_model_omega, breaks = 10000)
    hist(aln.data.LRT$alternative_model_omega, xlim = c(0,25), breaks = 10000)
    
    hist(aln.data.LRT$LRT_P_value_full) # number of alignments with significant LRT
    aln.data.LRT$LRT_P_value_log10[which(aln.data.LRT$LRT_P_value_log10 == -Inf)]
    table(cut(aln.data.LRT$LRT_P_value_log10, breaks = c(-Inf, seq(-9, 0, by = 1)), include.lowest = T, right = T))
    print(cumsum(table(cut(aln.data.LRT$LRT_P_value_log10, breaks = c(-Inf, seq(-9, 0, by = 1)), include.lowest = T, right = T))))
    
    sapply(names(aln.data.LRT)[c(5,7,8)], function(x) {
        print(length(which(aln.data.LRT[,x] < 0.01)));
    })
    sapply(names(aln.data.LRT)[c(5,7,8)], function(x) {
        print(length(which(aln.data.LRT[,x] < 0.05)));
    })
    
    # aln.p.value.cutoff <- 0.05
    # aln.data.LRT.significant.ensembl.id <- aln.data.LRT[which(aln.data.LRT$LRT_P_value_full < aln.p.value.cutoff),"ensembl.id"]
    # write.table(file = paste(c(file.alignment_codeml_results, ".LRT_5percsignificant_ensembl_id"), sep = "", collapse = ""),
                # aln.data.LRT.significant.ensembl.id, quote = F, sep = "\t", row.names = F, col.names = F)
    
    
    ### WRITE LRT RESULTS TABLE
    LTR.results.outfile <- paste("codeml_results_combined",
                                  "/",
                                 codeml.parameter.combination,
                                 ".alignment_codeml_results",
                                  ".LRT_results",
                                 sep = "")
                                               
    # write.table(file = LTR.results.outfile, aln.data.LRT, quote = F, sep = "\t", row.names = F, col.names = T)
    
    return(aln.data.LRT)
}

####### RESIDUE DATA #######
analyze.residue_codeml_results <- function(codeml.parameter.combination){

    ### read data
    file.residue_codeml_results <- paste("codeml_results_combined",
                                           "/",
                                           "codeml_results_base_files",
                                           "/",
                                           codeml.parameter.combination,
                                           ".residues_codeml_results",
                                           sep = "")
    print(file.residue_codeml_results)
    res.data <- read.table(file = file.residue_codeml_results, header = F, sep = "\t")
    
    names(res.data) <- c("ensembl.id",
                         "aa.human",
                         "position.aln",
                         "position.human",
                         "P.omega.gt.1",
                         "P.omega.gt.1.significance",
                         "omega.gt.1")
    head(res.data)
    dim(res.data)
    str(res.data)
    print(summary(res.data))
    

    ### analysis ###
    # aggregate(res.data, by = list(res.data$ensembl.id), FUN = length)
    res.data.ensembl.id.occurences <- ddply(res.data, "ensembl.id", summarise, N = length(ensembl.id))
    hist(log10(res.data.ensembl.id.occurences$N), breaks = 100)
    print(table(res.data.ensembl.id.occurences$N))
    
    barplot(table(res.data$aa.human))
    
    hist(res.data$P.omega.gt.1)
    print(table(res.data$P.omega.gt.1.significance))
    
    print(table(res.data[which(res.data$P.omega.gt.1 > 0.95),"P.omega.gt.1.significance"])) # ensembl.id with residues > 0.95 Pposterior
    print(table(res.data[which(res.data$P.omega.gt.1 > 0.99),"P.omega.gt.1.significance"]))
    
    print(nrow(res.data[which(res.data$P.omega.gt.1 > 0.95),]))
    print(length(sort(unique(res.data[which(res.data$P.omega.gt.1 > 0.95),"ensembl.id"]))))
    res.data.Ppos.significant.ensembl.id <- sort(unique(res.data[which(res.data$P.omega.gt.1 > 0.95),"ensembl.id"]))
    
#     write.table(file = paste(c(file.residue_codeml_results, ".Pposterior_95perc_ensembl_id"),
#                              sep = "", collapse = ""),
#                 res.data.Ppos.significant.ensembl.id, quote = F, sep = "\t", row.names = F, col.names = F)

    return(res.data)
}

####### COMBINED ALIGNMENT AND RESIDUE DATA #######
analyze.LRT_and_BEB_combined <- function(codeml.parameter.combination, aln.data.LRT, res.data){

    ### filter BEB results for 1) alignments that meet LRT < 0.05; and 2) Pposterior > 0.95
    
    data.LRT_and_BEB.merged <- merge(aln.data.LRT, res.data, by = "ensembl.id")
    dim(data.LRT_and_BEB.merged)
    head(data.LRT_and_BEB.merged)
    str(data.LRT_and_BEB.merged)
    
    LRT.P.cutoff <- 0.05
    BEB.P.cutoff <- 0.99
    
    data.LRT_and_BEB.merged.significant <-
        data.LRT_and_BEB.merged[which(data.LRT_and_BEB.merged$LRT_P_value_bh < LRT.P.cutoff &
                                      data.LRT_and_BEB.merged$P.omega.gt.1 > BEB.P.cutoff),]
    head(data.LRT_and_BEB.merged.significant)
    print(length(unique(data.LRT_and_BEB.merged.significant$ensembl.id)))
    
    
    ### write table
    LRT_and_BEB.results.outfile <- paste("codeml_results_combined",
                                 "/",
                                 codeml.parameter.combination,
                                 ".LRT_and_BEB_significant",
                                 sep = "")
    
    write.table(file = LRT_and_BEB.results.outfile, data.LRT_and_BEB.merged.significant, quote = F, sep = "\t", row.names = F, col.names = T)
    
    return(data.LRT_and_BEB.merged.significant)
}
###########################################################################


###########################
####### MAIN SCRIPT #######
###########################
codeml.parameter.combinations <- c("M7vM8_F61",
                                   "M7vM8_F3X4",
                                   "M1avM2a_F61",
                                   "M1avM2a_F3X4")
data.LRT_and_BEB.merged.significant.list <- list()
merge.columns <- c("ensembl.id", "aa.human", "position.aln", "length_alignment", "position.human", "length_human_sequence")

### for each of the codeml parameter combinations, read in the alignment-level and residue-level results and process them
for(codeml.parameter.combination in codeml.parameter.combinations){
    print(codeml.parameter.combination)
    
    aln.data.LRT <- analyze.alignment_codeml_results(codeml.parameter.combination)
    res.data <- analyze.residue_codeml_results(codeml.parameter.combination)
    data.LRT_and_BEB.merged.significant <- analyze.LRT_and_BEB_combined(codeml.parameter.combination, aln.data.LRT, res.data)
    
    # rename columns for later merging
    names(data.LRT_and_BEB.merged.significant)[which(!names(data.LRT_and_BEB.merged.significant) %in% merge.columns)] <- 
        paste(codeml.parameter.combination, "__", names(data.LRT_and_BEB.merged.significant)[which(!names(data.LRT_and_BEB.merged.significant) %in% merge.columns)], sep = "")
    
    names(data.LRT_and_BEB.merged.significant)
    
    # reorder columns
    data.LRT_and_BEB.merged.significant <- data.LRT_and_BEB.merged.significant[,c(1,13,14,2,15,3,4:12,16:18)]
    head(data.LRT_and_BEB.merged.significant)
    
    # store into list
    data.LRT_and_BEB.merged.significant.list[[codeml.parameter.combination]] <- data.LRT_and_BEB.merged.significant
}

### merge results from all parameter combinations: keep only those residues that are reported in all four
lapply(data.LRT_and_BEB.merged.significant.list, head, 2)
lapply(data.LRT_and_BEB.merged.significant.list, nrow)

# merge
data.LRT_and_BEB.merged.significant.all_parameters_combinations <- join_all(data.LRT_and_BEB.merged.significant.list, by = merge.columns, type = "inner")

# reorder
data.LRT_and_BEB.merged.significant.all_parameters_combinations<- data.LRT_and_BEB.merged.significant.all_parameters_combinations[order(data.LRT_and_BEB.merged.significant.all_parameters_combinations$ensembl.id),]

dim(data.LRT_and_BEB.merged.significant.all_parameters_combinations)
head(data.LRT_and_BEB.merged.significant.all_parameters_combinations)

length(unique(data.LRT_and_BEB.merged.significant.all_parameters_combinations$ensembl.id))
nrow(unique(data.LRT_and_BEB.merged.significant.all_parameters_combinations[,1:4]))

### write table
combined.results.outfile <- paste("codeml_results_combined",
                             "/",
                             "codeml_results_combined__LRT_and_BEB_significant__residues.txt",
                             sep = "")

write.table(file = combined.results.outfile, data.LRT_and_BEB.merged.significant.all_parameters_combinations, quote = F, sep = "\t", row.names = F, col.names = T)
