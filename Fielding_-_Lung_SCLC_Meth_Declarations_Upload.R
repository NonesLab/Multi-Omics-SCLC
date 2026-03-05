# **********
# Declares the directories, file names, variable names, colours, probe types and sets out options to be used for the study analysis
# **********

# **********
# Setting up directory for ChAMP import and results directory
# **********
array_type = "EPIC"
first_run = FALSE # make TRUE if this is the first time importing the idat files
samp_gp_name = "Sample_Group"
analysis = "SCLC_82Samples"
sub_analysis = "" # "longVshort_chemoIO", "longVshort_other", "chemIOVother_long", "chemIOVother_short" "longVshort_survival"

# **********
# Set up nomalisation options and potential covariates for removal of any batch effects 
# Acceptable normalisation methods = c("PBC"), c("BMIQ"), c("QN"), c("BMIQ", "QN")
# **********
normalise = FALSE # set to TRUE if you do not wish to use results from previous trial and have not previously run normalisation.
set_seed = TRUE
# if set_seed is FALSE a random seed between 1 and 5000 will be generated otherwise state the seed required in the statement below
seed = ifelse(set_seed, 553, ceiling(runif(1, 1, 5000)))
seed_txt = paste0("_seed_", seed)
normalisation_methods = c("BMIQ", "QN")
possible_covariates = c("Batch_Num") # nominate the variables within the pd (csv file with idat files) that could require batch correction

# **********
# Setting up differential methylation analysis
# **********
# 'sd" for standard deviation and 'mad' for median absolute deviation
dev_type = "sd" 
# nominate the minimum standard deviation of probes to be reported in heatmaps 
cutoff = 0.1 # 0.2 

# **********
# Setting up cutoff for dichotomising the data
# **********
dich_cutoff = 0.3

# **********
# nominate the directory and prefix for the study
# **********
study_prefix = "Fielding_-_Lung_All_Seq"
working_dir = paste0(study_dir, "Methylation/WGS_SCLC_Analysis/")
results_dir = paste0(working_dir, "Results/", analysis, "/")
if(!dir.exists(results_dir)){dir.create(results_dir)}
# Data store results
RData_results_dir = paste0(results_dir, "Data/RData/")
if(!dir.exists(RData_results_dir)){dir.create(RData_results_dir)}
# Differential Statistics results
Stats_results_dir = paste0(results_dir, "Differential_stats/")
if(!dir.exists(Stats_results_dir)){dir.create(Stats_results_dir)}
# RNA results
RNA_results_dir = paste0(results_dir, "RNA_Seq/")
if(!dir.exists(RNA_results_dir)){dir.create(RNA_results_dir)}
# Gene info results
Gene_results_dir = paste0(results_dir, "Gene_info/")
if(!dir.exists(Gene_results_dir)){dir.create(Gene_results_dir)}

# **********
# Make sure the column name for the sample ID in the sample information sheet, is labelled "Sample_Name" as per ChAMP
# **********
methylation_directory = paste0(gsub("/WGS_SCLC_Analysis", "", working_dir), "idats/SCLC_idats")

# **********
# Set up the directory where all Illumina probe reference data can be located
# **********
probe_dir = paste0(study_dir, "Methylation/Reference_Material/") 
if (array_type == "EPIC"){
  Probe_file= paste0(probe_dir, "EPIC_data/infinium-methylationepic-v-1-0-b5-manifest-file.csv")
  Probe_info = fread(Probe_file, skip=7)
  Probe_info = as.data.frame(Probe_info)
} else {
  Probe_file = paste0(probe_dir, "Meth450/HumanMethylation450_15017482_v1-2.csv")
  Probe_info = fread(Probe_file, skip=7)
  Probe_info = as.data.frame(Probe_info)
}
centromere_info = read.csv(paste0(probe_dir, "CHR_Centromere_Location.csv"), header=TRUE, check.names=FALSE)

# **********
# Set up ChAMP EPIC data directories (Zhou et al. & Nordlund et al.) for extra filtering
# **********
# Zhou 2017 manifest
EPIC_hg19_file = paste0(probe_dir, "Zhou_2017/EPIC.hg19.manifest.tsv")
# McCartney et al. 2016 SNP, rs_ch & multi mapper data
SNP_list = paste0(probe_dir, "McCartney_2016/SNP_EPIC_Array_List.txt")
rs_ch_list = paste0(probe_dir, "McCartney_2016/rs_ch_probes.csv")
Multimapper_list = paste0(probe_dir, "McCartney_2016/Multi_mappers_EPIC_array.txt")

# **********
# nominate the the study name for plot titles
# **********
study_name = "WGS SCLC DNA Methylation"

# **********
# nominate the different categories and levels to be reported
# **********
sample_sheet_cat_all = c("Sample_Group", "Batch_Num", "Survival_Group", "Treatment", "Meth_Cluster") 

sample_categories = c("Sample_Group", "Batch_Num", "Survival_Group", "Treatment", "Meth_Cluster") 

sample_levels = list(c("Adeno", "NSCLC", "Squamous", "sclc", "LC_neuro", "other_met", "benign", "Diff-quik"), 
                     c("Batch_1", "Batch_2", "Batch_3", "Batch_4", "Batch_5", "Batch_6", "Batch_7", "Batch_8", "Batch_9", "Batch_10"),
                     c("Long", "Short", "Other"),
                     c("chemo/immuno", "other"),
                     c("Group1", "Group2", "Group3", "Group4")) 

names(sample_levels) = sample_categories
names(sample_categories) = sample_categories

# **********
# nominate the colours to use to represent the different categories to be reported
# ensure each level per category is represented (i.e. number of colours should equal the total number of levels across all categories)
# **********
na_colour = "#FAFAFA"
na_border = "#AAAAAA"
sample_colours = list(# (Sample_Group)
                      c("#f7c0a6", "#e5fab7", "#c3d581", "khaki2", "#aad4da", "#f2beba", "#bcafd0", "#8dd38f"),
                      # (Batch_Num)
                      c(brewer.pal(6, "Pastel2"), brewer.pal(7, "Pastel1")[c(1,4,7)], "#cccccc"),  
                      # (long, short, other)
                      c("#8dd38f", "#f2beba", "khaki2"),
                      # chemo/immuno, other
                      c("orange", "dodgerblue"),
                      # Methylation clusters from RNASeq gene analysis of 156 genes
                      c("#b9e4c1", "#d3c3ef", "#f8d9b0", "#edeac4")
                      )                                      
sample_border_colours = list(# (Sample_Group)
                             c("coral", "darkolivegreen1", "darkolivegreen4", "goldenrod1", "cadetblue", "indianred1", "#68508b", "#349336"), 
                             # (Batch_Num)
                             c(brewer.pal(6, "Set2"), brewer.pal(7, "Set1")[c(1,4,7)], "#888888"),
                             # (long, short, other)
                             c("#349336", "indianred1", "goldenrod1"),
                             # chemo/immuno, other
                             c("darkorange2", "dodgerblue3"),
                             # Methylation clusters from RNASeq gene analysis of 156 genes
                             c("#49b05c", "#5108b4","#f58919", "#d9c80d")
                             )
names(sample_colours) = names(sample_border_colours) = sample_categories   
for (cat in names(sample_colours)){
  names(sample_colours[[cat]])=sample_levels[[cat]]
  names(sample_border_colours[[cat]])=sample_levels[[cat]]
}

# **********
# nominate the max number of probes to be reported in heatmaps and QC
# **********
numPos = 20000

# *************
# CHECK locations and filenames for the summary files to be read in
# *************
analysis_file = paste0(study_dir, "Extra_data/", study_prefix, "_Study_Report_GRCh37_June2024.csv")
WGS_VEP_data = paste0(study_dir, "GeneLandscape_Summary/WGS_SCLC_Analysis/RData/Fielding_-_Lung_All_Seq_VEP_76_WGS_filtered.RData")
# *************
# Reads in the nominated summary files
# *************
data.sets = read.in.data(analysis.file=analysis_file)
study_data = data.sets$study
study_data = create.unique.label(ending.type="all",
                                 study.data=study_data)
study_data %<>% arrange(UniqueLabel)
# remove duplicated records
study_data = study_data[which(duplicated(study_data$analysis)==FALSE),]
# remove diff quik, test and contaminated samples from the cohort
study_data %<>% filter(grepl("D0", donorLabel),
                       !testSampleLabel=="3300067",
                       !grepl("diluted", testSampleLabel),
                       !grepl("DQ", UniqueLabel))


