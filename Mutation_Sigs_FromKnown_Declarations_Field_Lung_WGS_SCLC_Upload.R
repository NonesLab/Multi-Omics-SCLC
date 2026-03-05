# ************************************************************** 
# Declaration sheet for Mutation Signatures from known source
# this sheet is to be edited where applicable 
# ************************************************************** 

# *************
# EDIT which signatures for detailed reporting 
# (must be one of either 'sample' (sample), 'consensus'(consensus call between sample & cohort) or 'cohort' (whole cohort))
# *************
report = "sample"

# *************
# EDIT the sample identifier in the study report to be displayed in plot and tables.
# *************
samp_id = "donorLabel"

# *************
# EDIT the study name to be used in plot and table titles
# *************
study_name = "Fielding Lung - SCLC WGS"

# *************
# EDIT the file prefix for the summarised study files and share directory
# ************* 
study_prefix = "Fielding_-_Lung_All_Seq"
working_dir = paste0(study_dir, "Signatures/WGS_SCLC_Analysis/")
results_dir = paste(working_dir, "Results/YAPSA/", sep="")

# *************
# EDIT what part of the study is to be summarised
# ************* 
# TRUE only if all samples (WGS & WES) in study are to be included
full_study = FALSE 
# if full_study is FALSE indicate whether WGS or WES is to be summarised
WES_only = FALSE
WGS_only = TRUE

# *************
# EDIT a temporary directory where temporary files can be written during the analysis
# *************
temp_dir = paste(working_dir, "Temp_Data/", sep="")

# *************
# CHECK location and filename for the study summary files to be read in
# *************
analysis_file = paste(study_dir, "Extra_data/", study_prefix, "_Study_Report_GRCh37_June2024.csv", sep="")

# *************
# CHECK location and filename for the study classifications to be read in
# *************
report_subgroups = TRUE
# If the above is true fill in the following variables where subgroup information can be found
sample_data_dir = paste0(study_dir, "Extra_data/WGS_SCLC_Analysis/")
sample_data_filename = "Fielding_-_Lung_SCLC_Meth_WGS_Samples.csv"
subgp_colname = "Meth_Cluster" 
samp_id = "donorLabel"

# *************
# read in the nominated analysis file and create sample data
# *************
all_data = read.in.study.data(pair.read=analysis_file,
                              full.study=full_study,
                              WGS.study=WGS_only,
                              WES.study=WES_only,
                              report = report_subgroups,
                              samp.dir.name=sample_data_dir,
                              samp.filename=sample_data_filename,
                              subgp.var=subgp_colname,
                              samp.id = samp_id)
study_data = all_data$study
sample_data = all_data$sample
# ***************
# CHECK the file ending for for the vcf files for samples
# ***************
# file endings for files where they exists in the somatic dna summary study report analysis filepath
vcf_file_ending = "/*.sp.vcf.gz"   

# *************
# EDIT location and parameters for signature matrix to be used for the analysis
# *************
# *************
# Choose signature matrix and cutoff for analysis 
# Two matrices from YAPSA package have been included here
# Either choose one of COSMIC v2 30 or COSMIC v3 47 
# OR comment out code below and supply your own matrix and assign to parameter 'sig_matrix' with colnames in form 'Sig.n'
# *************
# signature_matrix_name = "COSMIC v2 (30 Sig)"
# signature_matrix_name = "COSMIC v3 Real"
# signature_matrix_name = "COSMIC v3 Artif"
signature_matrix_name = "COSMIC v3.2 SBS"

# indel_sig_matrix_name = "PCAWG Indel v3 (16 Sig)"
indel_sig_matrix_name = "COSMIC v3.2 Indels"
sig_cutoff = 0.06 
indel_sig_cutoff = 0.06 

# *************
# EDIT reference genome to be used for YAPSA package
# *************
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
ref_genome = BSgenome.Hsapiens.UCSC.hg19


#######################################################################################################
# -----------This section is to set up all the script variables from the above set input---------------
#######################################################################################################
# *************
# create directories if they don't exisst
# *************
create.dir(results_dir)
create.dir(temp_dir)

# *************
# load signature details for the analyses
# *************
sig_details = get.sig.matrix.details(sbs.matrix.name = signature_matrix_name,
                                     indel.matrix.name = indel_sig_matrix_name,
                                     sbs.cutoff = sig_cutoff,
                                     indel.cutoff = indel_sig_cutoff,
                                     package="YAPSA")
# Assign details to variables
sig_matrix = sig_details$sbs_matrix
sig_cutoff_det = sig_details$cutoff
sig_ind_df = sig_details$sbs_sig_df
sig_order = sig_details$sbs_order
sig_colours = sig_details$sbs_colour
indel_sig_matrix = sig_details$indel_matrix
indel_sig_cutoff = sig_details$indel_cutoff
indel_sig_ind_df = sig_details$indel_sig_df
indel_colours = sig_details$indel_colour
if (!"category" %in% colnames(sig_ind_df)){
  sig_ind_df$category = ""
  for (cat in names(sig_details$combined_cols)){
    if (cat =="Age"){
      sig_ind_df$process[which(sig_ind_df$sig=="Sig.5")] = "clock-like signature"
      sig_ind_df$category[which(sig_ind_df$sig%in%c("Sig.1", "Sig.5"))] = cat
    } else if (cat=="HR_Def") {
      sig_ind_df$category[grep("HR|repair ", sig_ind_df$process)] = cat
    } else if (cat=="MMR") {
      sig_ind_df$category[intersect(grep("MMR", sig_ind_df$process), grep("POL", sig_ind_df$process, invert=TRUE))] = cat
    } else if (cat=="POLE") {
      sig_ind_df$category[intersect(grep("POL E|Eta", sig_ind_df$process, ignore.case = TRUE), grep("MMR", sig_ind_df$process, invert=TRUE))] = cat
    } else if (cat=="POLD1") {
      sig_ind_df$category[intersect(grep("POLD1", sig_ind_df$process), grep("MMR", sig_ind_df$process, invert=TRUE))] = cat
    } else if (cat=="POLE_MMR") {
      sig_ind_df$category[intersect(grep("POL E|Eta", sig_ind_df$process, ignore.case = TRUE), grep("MMR", sig_ind_df$process))] = cat
    } else if (cat=="POLD1_MMR") {
      sig_ind_df$category[intersect(grep("POLD1", sig_ind_df$process), grep("MMR", sig_ind_df$process))] = cat
    } else if (length(grep("ROS", cat))>0) {
      sig_ind_df$category[grep("Reactive oxygen", sig_ind_df$process)] = cat
    } else if (cat=="Haloalkane") {
      sig_ind_df$category[grep("Halogenalkane", sig_ind_df$process)] = cat
    } else if (cat=="Artifact") {
      sig_ind_df$category[grep("nonvalidated", sig_ind_df$process)] = cat
    } else if (length(grep("Aristo", cat))>0) {
      sig_ind_df$category[grep("Aristo", sig_ind_df$process, ignore.case = TRUE)] = "Aristolochic"
    } else {
      sig_ind_df$category[grep(cat, sig_ind_df$process, ignore.case = TRUE)] = cat
    }
  }
  sig_ind_df$category[which(sig_ind_df$category=="")] = "Unknown"
}
sig_ind_df %<>% mutate(cat_cols=sapply(category, function(x) ifelse(x=="ROS",
                                                                    "#225EA8",
                                                                    sig_details$combined_cols[x])))
if (!"category" %in% colnames(indel_sig_ind_df)){
  indel_sig_ind_df$category = ""
  for (cat in names(sig_details$combined_cols)){
    if (cat =="Age"){
      indel_sig_ind_df$category[which(indel_sig_ind_df$sig%in%c("ID1", "ID2"))] = cat #paste(cat, "MMR", sep="_")
    } else if (cat=="HR_Def") {
      indel_sig_ind_df$category[grep("Defective homologous recombincation DNA damage repair|defective HR repair", indel_sig_ind_df$process)] = cat
    } else if (cat=="MMR") {
      indel_sig_ind_df$category[grep("Defective DNA mismatch repair", indel_sig_ind_df$process)] = cat
    } else if (cat=="UV") {
      indel_sig_ind_df$category[grep("Ultraviolet light", indel_sig_ind_df$process)] = cat
    } else if (cat=="POLE_MMR") {
      if (length(grep("Repair of DNA double strand breaks by NHEJ", indel_sig_ind_df$process))>0){
        indel_sig_ind_df$category[grep("Repair of DNA double strand breaks by NHEJ", indel_sig_ind_df$process)] = "NHEJ_TOP2A"
      } else {
        indel_sig_ind_df$category[intersect(grep("DBS repair by non-homologous end join", indel_sig_ind_df$process),
                                            grep("defective HR repair", indel_sig_ind_df$process, invert = TRUE))] = "NHEJ"
      }
    } else if (cat=="POLD1_MMR") {
      indel_sig_ind_df$category[grep("topoisomerase TOP2A", indel_sig_ind_df$process)] = "TOP2A"
    } else {
      indel_sig_ind_df$category[grep(cat, indel_sig_ind_df$process, ignore.case = TRUE)] = cat
    }
  }
}
indel_sig_ind_df %<>% mutate(cat_cols=sapply(category, function(x) ifelse(length(grep("NHEJ", x))>0,
                                                                          sig_details$combined_cols["POLE_MMR"],
                                                                          ifelse(x=="TOP2A",
                                                                                 sig_details$combined_cols["POLD1_MMR"],
                                                                                 ifelse(x=="Age", # "Age_MMR",
                                                                                        sig_details$combined_cols["Age"],
                                                                                        sig_details$combined_cols[x])))))
if (subgp_colname=="Meth_Cluster"){
  level_cols = level_cols = c("#49b05c", "#5108b4","#f58919", "#d9c80d", "#f188d5", "#610aa3", "#e62843")
} else if (subgp_colname=="Cat"){
  level_cols = c("coral", "darkolivegreen1", "darkolivegreen4", "lightgoldenrod1", "cadetblue", "indianred1", "mediumorchid3", "darkseagreen1", "#555555") 
} 

# ***********************
# for analysis name the signature packages to use
# ***********************
sig_packages = c("YAPSA") 
sig_cutoff = list()
sig_cutoff_text = c()
for (package in sig_packages){
  if (package=="YAPSA"){
    sig_cutoff[["YAPSA"]] = sig_cutoff_det
    sig_cutoff_text = c(sig_cutoff_text, ifelse(length(sig_cutoff$YAPSA)==1, sig_cutoff$YAPSA, "optimised"))
  } 
}
if (WES_only) {
  if (length(indel_sig_cutoff)>1){
    indel_sig_cutoff = 0.1
  }
}
