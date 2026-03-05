# ************************************************************** 
# Declaration sheet for GeneLandscape Reports
# i.e Singe Gene, Nominated Gene list or extracting top mutated genes for study
# this sheet is to be editted where applicable 
# ************************************************************** 
# ***************
# EDIT the directory and prefix for the study
# ***************  
study_prefix = "Fielding_-_Lung_All_Seq" 
working_dir = paste0(study_dir, "GeneLandscape_Summary/")

# ***************
# EDIT the analysis and study names
# ***************
study_name = "Fielding Lung - SCLC WGS" # "Fielding Lung Illumina Sequencing" "Fielding Lung WES" "TCGA Adenocarcinoma WES" 
analysis_name = "Gene Variants Landscape"

# *************
# CHECK location and filename for the study summary files to be read in
# *************
analysis_file = paste(study_dir, "Extra_data/", study_prefix, "_Study_Report_GRCh37_June2024.csv", sep="")

# *************
# EDIT which part of the study is to be summarised
# ************* 
full_study = FALSE # mark as true if both WGS and WES are to be summarised for study
# if not a full study you must mark one of WGS or WES as TRUE
WGS_only = TRUE
WES_only = FALSE

# *************
# EDIT if list of genes is to be generated or genes are to be summarised
# ************* 
first_pass = FALSE # mark as TRUE if a gene list is to be generated and copied to results directory

# *************
# EDIT the reporting that will be used
# ************* 
report_just_loss = TRUE # TRUE = del & cnLOH will be reported as Loss in summary matrix
report_just_HG = TRUE   # TRUE = only HighGains (i.e. >=CN9 for ploidy >=2.7 & >=C5 for ploidy <2.7) are reported in summary matrix

# *************
# Annotations for reporting in heatmaps
# ************* 
# MUST add Status (i.e. tumour type) to list where you want it reported
top_annots = c("Amp_Arch", "WGD_Status", "MSI", "TMB", "Ploidy", "Treatment", "Survival", "Cellularity", "Meth_Cluster", "Status")  
# mark as TRUE if cases reported are to be from the DonorLabel otherwise FALSE will report unique label, UniqSampleLabel
by_donor = TRUE
# mark as TRUE if reporting sample status/subgroups in plots and nominate file below with information
status = TRUE 
# mark as TRUE if a specific sample order is to followed for all plots
sample_order = TRUE

# If the above is true fill in the following variables where subgroup information can be found
sample_data_dir = paste0(study_dir, "Extra_data/WGS_SCLC_Analysis/")
sample_data_filename = "Fielding_-_Lung_SCLC_Meth_WGS_Samples_v2.csv"
subgp_colname = "Cat"
samp_id = "donorLabel"

# ***************
# EDIT if first pass is FALSE and/or a list of genes for the study has been established
# ***************
# Nominate report type
report_type = "GeneList" # "Single"  
report_gene = "Drivers" 
# ******************************
# if wanting to look at common ctDNA mutations or copy number/LOH then nominate the genes to be reported otherwise set to NULL (a list of genes must still be defined)
top_genes = NULL # c("TP53", "RB1") # 
# mark as TRUE to compare with ctDNA and FALSE to report copy number/LOH
add_ctdna = FALSE
# ******************************
if (!first_pass){
  #### Gene List ####
  if (report_type=="GeneList"){
    # suppy location of gene list NOTE: Must be identified with Hugo_Symbol or single gene.
    if (report_gene == "Drivers"){
      gene_file = paste(working_dir, "Spreadsheets/Genes_To_Report/Fielding_-_Lung_George_SCLC_genes_withGrouping.csv", sep="") ############# _withGrouping ##############
    } 
  } else {
    gene_file = ""
  }
  #### Single Gene ####
  if (report_type=="Single"){
    # if single gene study nominate the gene to be summarised 
    Gene = "EGFR"
  } else {
    Gene = ""
  }
}

# ***************
# EDIT if first pass is FALSE when running full summary report on study
# ***************
# ***************
# Nominate flags for the type of report to be generated for the study
# ***************
specific_gene_order = FALSE
if (specific_gene_order){
  if (report_gene=="Drivers"){
    full_gene_order_file = ""
    SNV_gene_order_file = NULL
  } 
} else {
  full_gene_order_file = NULL
  SNV_gene_order_file = NULL
}
#  If status is to be included set up the status file
create.dir(working_dir)
if (status){
  Status_data = as.data.frame(fread(paste(sample_data_dir, sample_data_filename, sep=""))) 
  if (length(grep("UniqueLabel", colnames(Status_data)))==0){
    Status_data %<>% mutate(UniqueLabel = paste(Status_data$donorLabel, "WGS", sep="_"))
  }
  if (sample_order){
    sample_levels = dplyr::arrange(Status_data, Meth_Index)
  }
} else {
  Status_data = NULL
  if (sample_order){
    sample_levels = as.data.frame(fread(paste(gsub("GeneLandscape_Summary", "Extra_data", working_dir), "WGS_SCLC_Analysis/", sample_data_filename, sep="")))
    sample_levels %<>% dplyr::arrange(Meth_Index)
  } else {
    sample_levels = NULL
  }
}

#######################################################################################################
# -----------This section is to set up all the script variables from the above set input---------------
#######################################################################################################
# *************
# create directories if they don't exist
# *************
results_dir = paste(working_dir, "WGS_SCLC_Analysis/", study_prefix, ifelse(report_type=="Single", paste("_", Gene, sep=""), ""), "_MethSamp_Results/", sep="")
temp_dir = paste(working_dir, "temp/", sep="")
create.dir(results_dir)
create.dir(temp_dir)

# *************
# Reads in the nominated summary files and create the subgroups for the analyses
# *************
data.sets = read.in.study.data(analysis.file=analysis_file,
                               full.study=full_study,
                               WGS.study=WGS_only,
                               WES.study=WES_only,
                               report = status,
                               samp.dir.name=sample_data_dir,
                               samp.filename=sample_data_filename,
                               subgp.var=subgp_colname,
                               samp.id = samp_id,
                               end.type="all")
study_data = data.sets$study
sample_data = data.sets$sample
Status_data = data.sets$status
# for SCLC WGS project
Status_data %<>% dplyr::mutate(WGS_Nonsyn_TMB = unlist(sapply(UniqueLabel, function(x) study_data$TMB_nonSyn[which(study_data$UniqueLabel==x)])))
#  Add all the annotations to the study data
case_to_report = "UniqueLabel"
annot_data = add.annotations(study.data=study_data,
                             top.annots = top_annots,
                             study.prefix = study_prefix,
                             work.dir = working_dir,
                             case = case_to_report,
                             status = status,
                             status.file = Status_data,
                             subgp.col = subgp_colname)
study_data = annot_data$study
study_data[,samp_id] = factor(study_data[,samp_id], levels = sample_levels[,samp_id])
study_data = study_data[order(study_data[,samp_id]),]

# ***************
# Enter mutation types to be reported from list of possible variants across the study (i.e. first_run = FALSE)
# and enter the number of genes you wish reported on (i.e. first_pass = TRUE)
# ***************
if (!first_pass) {
  # create donor list
  if (sample_order){
    Donor_list = sample_levels[,samp_id]
  } else {
    Donor_list = as.character(unique(study_data[,samp_id]))
  }
  # ***************
  # establish the gene list and the donor list for the rest of the analysis
  # ***************
  if (report_type=="GeneList"){
    Genes_of_interest = as.data.frame(fread(gene_file, header=TRUE, sep=",", na.strings="",
                                            check.names=FALSE, stringsAsFactors=FALSE))
    if ("Group" %in% colnames(Genes_of_interest)){
      # Need to add ZFHX3 to the list 
      Genes_of_interest = rbind(Genes_of_interest,
                                data.frame(Hugo_Symbol="ZFHX3",
                                           SNV = 2,
                                           CNV = 59,
                                           SV=2,
                                           Other_Name = NA,
                                           Group = "Transcriptional_Regulation"))
      Gene_list = data.frame(gene_num = c(1:nrow(Genes_of_interest)),
                             Hugo_Symbol = Genes_of_interest$Hugo_Symbol,
                             stringsAsFactors = FALSE)
      Genes_of_interest$Group = factor(Genes_of_interest$Group, levels = c("Cell_Cycle_Regulation", "Transcriptional_Regulation", "Receptor_Kinase_P13K_signalling", "Notch_signalling_NE_differentiation", "Other"))
      Genes_of_interest %<>% dplyr::arrange(Group, desc(SNV))
      Gene_list = data.frame(Gene_list,
                             Group = Genes_of_interest$Group)
      Gene_list$Group = factor(Gene_list$Group, levels = c("Cell_Cycle_Regulation", "Transcriptional_Regulation", "Receptor_Kinase_P13K_signalling", "Notch_signalling_NE_differentiation", "Other"))
    } else {
      Genes_of_interest %<>% dplyr::select(Hugo_Symbol, SNV, CNV, SV) %>% dplyr::filter(Hugo_Symbol%in%c("TP53", "RB1", "EP300", "CREBBP", "FMN2", "TP73"))
      Gene_list = dplyr::mutate(Genes_of_interest, gene_num = c(1:6)) %>% dplyr::select(gene_num, Hugo_Symbol)
    }
    Gene_list$Hugo_Symbol = factor(Gene_list$Hugo_Symbol, levels=c(as.character(Gene_list$Hugo_Symbol)))
  }
} 

# ***************
# list of the existing SNVs to be reported in the summary
# ***************
VEP_anno = TRUE
if (VEP_anno){
  SNV_mutations = c("transcription_ablation", "splice_acceptor", "splice_donor", "splice_region", "stop_gained", "start_gained", 
                    "frameshift_INS", "frameshift_DEL", "stop_lost", "start_lost", "transcription_amplification", "inframe_INS", "inframe_DEL",
                    "missense", "protein_altering", "incomplete_terminal_codon", "regulatory_region_ablation", "germline") # "stop_retained", "start_retained", 
} else {
  SNV_mutations = c("Nonsense_Mutation", "Frame_Shift_INS", "Frame_Shift_DEL", "In_Frame_Ins", "In_Frame_Del",
                    "Missense_Mutation", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation", "germline")
}
# ***************
# list of the existing CNs to be reported in the summary
# ***************
CNV_mutations = c("Del", "cnLOH", "Gain", "HighGain", "Het")
CNV_reported = c("Loss", "TotalLoss", "cnLOH", "Gain 3-4", "Gain 5-8", "HighGain 5-8", "HighGain 9=<")

# ***************
# list of the SVs to be reported in the summary
# ***************
SV_mutations = c("Loss_of_Function", "Fusion")
SNV_file = "sp.maf.gz"  # (from analysis directory)
SV_WGS_file = "*somatic.ascat.mere.type.annotated.HighConfidence.dcc"
SV_WES_file = "*somatic.filtered.sequenza.mere.type.annotated.HighConfidence.dcc"
CNV_file = ".genes.pc.best.cn.txt"   # from either ascat or sequenza analysis directory
