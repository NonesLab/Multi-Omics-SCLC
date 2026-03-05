# ***************
# This script is a repository of all the functions
# required when running the GeneLandscape summary script for the Fielding Lung sequencing project
# ***************
# ***************
# Function to create a directory if it doesn't exist
# ***************
create.dir = function(directory.name){
  if(!dir.exists(directory.name)){
    dir.create(directory.name, recursive=TRUE)
  }
}

# ***************
# Function to  read in the study summary spreadsheets for the mutations analysis summary
# ***************
read.in.data = function(analysis.file = NULL,
                        full.study=TRUE, 
                        WGS.study=FALSE,
                        WES.study=FALSE){
  if (!is.null(analysis.file)){
    study.data = as.data.frame(fread(file=analysis.file)) 
  } else {
    study.data = NULL
  }
  if (!is.null(study.data)){
   if (full.study){
      study.data %<>% dplyr::filter(!is.na(analysis)) %>% dplyr::filter(!analysis=="")
    } else {
      study.data %<>% dplyr::filter(!is.na(analysis)) %>% dplyr::filter(!analysis=="")
      if (WGS.study){
        study.data %<>% dplyr::filter(is.na(testCaptureKit)|testCaptureKit=="") %>% dplyr::filter(is.na(controlCaptureKit)|controlCaptureKit=="")
      } else if (WES.study) {
        study.data %<>% dplyr::filter(!is.na(testCaptureKit)) %>% dplyr::filter(!is.na(controlCaptureKit))
      }
    }
  }
  data.sets = list(study=study.data)
  return(data.sets)
} 

# ***************
# Function to build the file and directory names where the study, SNVs, CNVs and SVs are located
# ***************
build.filename = function(UniqLab,
                          study.data, 
                          ending){
  pos = which(study.data$UniqueLabel==UniqLab)
  filename = paste(study.data$donorLabel[pos], ".", 
                   study.data$testSampleLabel[pos], ".",
                   unlist(strsplit(study.data$testAlignedReadGroupSetIRI[pos], ":"))[2], "_vs_",
                   study.data$donorLabel[pos], ".",
                   study.data$controlSampleLabel[pos], ".",
                   unlist(strsplit(study.data$controlAlignedReadGroupSetIRI[pos], ":"))[2], ".",
                   ending, sep="")
  if (study.data$donorLabel[pos]=="D03_19_027"){
    filename = gsub("D03_19_", "D01_19_", filename)
  }
  return(filename)
}

# ***************
# Function to create a unique label for every sample in the analysis
# ***************
create.unique.label = function(by.donor = TRUE,
                               ending.type = "seq",
                               study.data = study_data){
  # ***************
  # Function to create ending for a unique label for every sample in the analysis
  # ***************
  BGI.seq = c("BGI", "BGIseq500", "MGIseq2000")
  Illumina.seq = c("Illumina", "NovaSeq", "HiSeq2000", "HiSeq2500", "HiSeqXTen", "MiSeq", 
                   "NextSeq", "HiSeqXFive", "HiSeq4000", "HiSeq1500", "HiSeq3000")
  get.ending = function(study.data=study.data,
                        sample.vector,
                        ending.type){
    pos = ifelse(ending.type == "seq" | ending.type == "all",
                 which(colnames(study.data)=="testCaptureKit"),
                 ifelse(ending.type == "platform",
                        which(colnames(study.data)=="testSequencingPlatform"),
                        ifelse(ending.type=="DQ_WGS",
                               which(colnames(study.data)=="studyAltLabel"),
                               which(colnames(study.data)=="donorLabel"))))
    # message(paste("The position of comparison is:", pos))
    ending = ifelse(ending.type == "seq" | ending.type == "all", 
                    ifelse(length(grep("custom", sample.vector[pos], ignore.case=TRUE))>0,
                           "Custom",
                           ifelse(length(grep("methyl", sample.vector[pos], ignore.case=TRUE))>0,
                                  "Methyl",
                                  ifelse(length(grep("panel|amplicon", sample.vector[pos], ignore.case=TRUE))>0,
                                         "Panel",
                                         ifelse(length(grep("exo", sample.vector[pos], ignore.case=TRUE))>0,
                                         "WES", "WGS")))),
                    ifelse(ending.type=="platform",
                           ifelse(as.character(sample.vector[pos])%in%Illumina.seq, 
                                  "Illumina",
                                  ifelse(as.character(sample.vector[pos])%in%BGI.seq, "BGI")), 
                           ifelse(ending.type=="DQ_WGS", 
                                  ifelse(as.character(sample.vector[pos])=="Fresh", "Fresh", "Diff"),"wrong")))
    # message(paste("The ending given to the comparison is:", ending))
    return(ending)
  }
  if (by.donor){
    if (unique(is.na(study.data$donorLabel)==TRUE)){
      by.donor = FALSE
      } else if (length(which(is.na(study.data$donorLabel)))>0){
        by.donor = FALSE
      }
  }
  if (by.donor){
    study.data %<>% mutate(UniqueLabel = apply(study.data, 1, function(x) paste(x[which(colnames(study.data)=="donorLabel")], 
                                                                                "_", get.ending(study.data, x, ending.type), sep="")))
  } else {
    study.data %<>% mutate(UniqueLabel = apply(study.data, 1, function(x) paste(x[which(colnames(study.data)=="testSampleLabel")], 
                                                                                "_", get.ending(study.data, x, ending.type), sep="")))
  }
  study.data$UniqueLabel[grep("dilute", study.data$testSampleLabel)] = paste(study.data$UniqueLabel[grep("dilute", study.data$testSampleLabel)], "D", sep="_")
  if (ending.type == "all"){
    for( i in 1:nrow(study.data)){
      study.data$UniqueLabel[i] = ifelse(study.data$testSequencingPlatform[i]%in%Illumina.seq, paste0(study.data$UniqueLabel[i], "_I"), paste0(study.data$UniqueLabel[i], "_B"))
      study.data$UniqueLabel[i] = ifelse(length(grep("3300", study.data$testSampleLabel[i]))>0,
                                         paste0(study.data$UniqueLabel[i], "FF"),
                                         ifelse(length(grep("5500", study.data$testSampleLabel[i]))>0,
                                                paste0(study.data$UniqueLabel[i], "DQ"),
                                                ifelse(length(grep("6600", study.data$testSampleLabel[i]))>0,
                                                       paste0(study.data$UniqueLabel[i], "CT"),
                                                       paste0(study.data$UniqueLabel[i], "FF"))))
    }
  }
  Labeldups = which(duplicated(study.data$UniqueLabel))
  if (length(Labeldups)>0){
    labels = unique(study.data$UniqueLabel[Labeldups])
    for (label in labels){
      dups = which(study.data$UniqueLabel==label)
      new_labels = paste(label, seq(1,length(dups),1), sep="_")
      for (i in 1:length(dups)){
        study.data$UniqueLabel[dups[i]] = new_labels[i]
      }
    }
  }
  return(study.data)
}

# *************
# Function to create subgroups for reporting
# *************
create.subgroup.data = function(report,
                                samp.dir.name,
                                samp.filename,
                                subgp.var,
                                study.data,
                                samp.id){
  sample.data = NULL
  if (report && !is.null(samp.dir.name) && !is.null(samp.filename)){
    sample.data = as.data.frame(fread(paste(samp.dir.name, samp.filename, sep=""))) 
  } else if (report && is.null(samp.dir.name) || is.null(samp.filename)) {
    stop("If you wish to report subgroups you must supply a directory and filename where the information can be found.")
  }
  if (!is.null(sample.data)){
    samp.col = names(which(apply(sample.data, 2, function(x) length(which(x %in% study.data[,samp.id]))>0)==TRUE))[1]
    subgroups = data.frame(subgroup = sample.data[, subgp.var],
                           PID = sample.data[,samp.col])
    subgroup.df = data.frame(PID = vector(mode="character"),
                             subgroup = vector(mode="character"))
    for (i in 1:nrow(subgroups)){
      sample.labels = study.data$UniqueLabel[grep(subgroups$PID[i], study.data$UniqueLabel)]
      for (label in sample.labels){
        subgroup.df = rbind(subgroup.df, data.frame(PID=label, 
                                                    subgroup=subgroups$subgroup[i]))
      }
    }
    if (length(which(is.na(subgroup.df$subgroup))>0)){
      message("All NAs found in the input subgroups were replaced with 'Other' ")
      subgroup.df$subgroup[which(is.na(subgroup.df$subgroup))] = "other"
    }
    subgroup.df$subgroup = factor(subgroup.df$subgroup, levels=c("adeno", "NSCLC", "squamous", "sclc", "LC_neuro", "other_met", "other"))
    level_cols = c("coral", "darkolivegreen1", "darkolivegreen4", "goldenrod1", 
                   "cadetblue3", "indianred1", "mediumorchid3", "darkseagreen1", "#555555") 
    subgroup.df$col = sapply(subgroup.df$subgroup, function(x) level_cols[which(levels(subgroup.df$subgroup)==x)])
  } else {
    subgroup.df = data.frame(PID = study.data$UniqueLabel,
                             subgroup = "Subgroup",
                             col = "#aaaaaa")
  }
  return(subgroup.df)
}


# *************
# Function to create the sample data with subgroups for reporting
# *************
create.sample.data = function(report,
                              samp.dir.name,
                              samp.filename,
                              subgp.var,
                              study.data,
                              study.bgi.data=NULL,
                              samp.id){
  samp.study.data = create.subgroup.data(report = report,
                                         samp.dir.name = samp.dir.name,
                                         samp.filename = samp.filename,
                                         subgp.var = subgp.var,
                                         study.data = study.data,
                                         samp.id = samp.id)
  return(samp.study.data)
}

# *************
# function to read in the nominated analysis file and create sample data
# *************
read.in.study.data = function(analysis.file=analysis_file,
                              full.study=full_study,
                              WGS.study=WGS_only,
                              WES.study=WES_only,
                              report = NULL,
                              samp.dir.name=sample_data_dir,
                              samp.filename=sample_data_filename,
                              subgp.var=subgp_colname,
                              samp.id = samp_id,
                              end.type="seq"){
  data.sets = read.in.data(analysis.file=analysis.file,
                           full.study=full.study, 
                           WGS.study=WGS.study,
                           WES.study=WES.study)
  study.data = data.sets$study
  study.data$donorLabel = gsub("-", "_", study.data$donorLabel)
  # *************
  # creating an unique sample identifier for reporting
  # *************
  study.data = create.unique.label(study.data=study.data,
                                   ending.type = end.type)
  # *************
  # Filter summary files to samples of interest
  # *************
  samp.of.interest = as.data.frame(fread(paste(samp.dir.name, samp.filename, sep="")))
  colnames(samp.of.interest)[grep("sequencing", colnames(samp.of.interest))] = "testSequencingPlatform"
  samp.of.interest = create.unique.label(study.data=samp.of.interest,
                                         ending.type = end.type)
  study.data %<>% filter(UniqueLabel%in%samp.of.interest$UniqueLabel) %>% filter(testSampleLabel%in%samp.of.interest$testSampleLabel)
  # remove diluted samples for 100 Illumina exome study
  pattern = paste(paste(paste("WES", seq(2,10,1), sep="_"), collapse="|"), paste(paste("WGS", seq(2,10,1), sep="_"), collapse="|"), sep="|")
  study.data %<>% filter(!grepl(pattern, UniqueLabel))
  study.data$UniqueLabel = gsub("WES_1", "WES", gsub("WGS_1", "WGS", study.data$UniqueLabel))
  study.data %<>% filter(!grepl("D01_19_004_WG", UniqueLabel))
  # *************
  # create the subgroups for the analyses
  # *************
  sample.data = create.sample.data(report = status,
                                   samp.dir.name=samp.dir.name,
                                   samp.filename=samp.filename,
                                   subgp.var=subgp.var,
                                   study.data=study.data,
                                   study.bgi.data=bgi.study.data,
                                   samp.id = samp.id)
  return(list(study=study.data, sample=sample.data, status=samp.of.interest))
}


# ***************
# Function to populate df with cohort SNV 
# ***************
populate.snv = function(study.data,
                        temp.dir=temp_dir){
  # ***************
  # Set up dataframe
  # ***************
  SNV.data = NULL
  SNV.stats = NULL
  # ***************
  # Populate dataframe
  # ***************
  for (i in 1:nrow(study.data)){
    sample.data = NULL
    # ***************
    # read in SNV data
    # ***************
    if (!is.na(study.data$analysisPath[i])){
      sample.file = list.files(path=study.data$analysisPath[i], pattern=build.filename(UniqLab = study.data$UniqueLabel[i],
                                                                                       study.data = study.data,
                                                                                       SNV_file), full.names=TRUE)
      # read in SNP data
      if (file.exists(sample.file)){
        message(paste(i, "SNV file found for", study.data$UniqueLabel[i]))
        del.header = as.integer(system(paste("zcat ", sample.file, " | grep -c '#' ", sep=""), intern=TRUE))
        system(paste("zcat ", sample.file, " | cut -f ", paste(seq(1, 62, 1), collapse=","), " > ", temp.dir, "temp.tsv", sep=""), intern=TRUE)
        sample.data = fread(paste(temp.dir, "temp.tsv", sep=""), skip=del.header,
                            na.strings="NA", strip.white=TRUE, sep="\t",
                            colClasses=c("character"), stringsAsFactors=FALSE)
        system(paste("rm ", temp.dir, "temp.tsv", sep=""))
      } else {
        message(paste("SNV file", sample.file, "NOT found for", study.data$UniqueLabel[i], "\n"))
      }
    }
    # ***************
    # create SNV_data
    # ***************
    if (!is.null(sample.data) && nrow(sample.data) != 0){
      sample.stats = data.frame(DonorLabel = study.data$donorLabel[i],
                                UniqueLabel = study.data$UniqueLabel[i],
                                Somatic_SNP = nrow(filter(sample.data, Variant_Type%in%c("SNP", "DNP", "TNP"))),
                                Somatic_dbSNP = nrow(filter(sample.data, Variant_Type%in%c("SNP", "DNP", "TNP"), !DbSNP_RS=="novel")),
                                Somatic_Indel = nrow(filter(sample.data, Variant_Type%in%c("INS", "DEL"))),
                                Somatic_dbIndel = nrow(filter(sample.data, Variant_Type%in%c("INS", "DEL"), !DbSNP_RS=="novel")))
      # message(paste("SNV data has been read in for", study.data$UniqueLabel[i]))
      sample.data %<>% mutate(DonorLabel = study.data$donorLabel[i],
                              UniqSampleLabel = study.data$UniqueLabel[i],
                              VAF = as.integer(T_Alt_Count)/as.integer(T_Depth)) %>%
        select(DonorLabel,
               UniqSampleLabel,
               Hugo_Symbol,
               Chromosome,
               Start_Position,
               End_Position,
               Strand,
               Reference_Allele,
               Tumor_Seq_Allele1,
               Tumor_Seq_Allele2,
               VAF,
               QFlag,
               Confidence,
               Mutation_Status,
               Variant_Classification,
               Variant_Type,
               Eff_Impact,
               Effect_Ontology,
               Effect_Class,
               Transcript_BioType,
               Transcript_ID,
               Amino_Acid_Change,
               CDS_Change,
               Tumor_Sample_Barcode)
      SNV.data = rbind(SNV.data, sample.data)
      SNV.stats = rbind(SNV.stats, sample.stats)
    }
  }
  if (VEP_anno){
    SNV.data$Plot_Class = gsub("_variant", "", gsub("disruptive_", "", SNV.data$Effect_Ontology))
    SNV.data$Plot_Class[which(SNV.data$Plot_Class=="5_prime_UTR_premature_start_codon_gain")] = "start_gained"
    SNV.data$Plot_Class[which(SNV.data$Plot_Class=="frameshift")] = paste(SNV.data$Plot_Class[which(SNV.data$Plot_Class=="frameshift")],
                                                                          SNV.data$Variant_Type[which(SNV.data$Plot_Class=="frameshift")], sep="_")
    SNV.data %<>% mutate(Plot_Class = apply(SNV.data, 1,
                                            function(x) ifelse(x[which(colnames(SNV.data)=="Plot_Class")]=="frameshift",
                                                               ifelse(x[which(colnames(SNV.data)=="Reference_Allele")]=="-",
                                                                      "frameshift_insertion",
                                                                      "frameshift_deletion"),
                                                               x[which(colnames(SNV.data)=="Plot_Class")])))
    SNV.data %<>% mutate(Plot_Class = apply(SNV.data, 1,
                                            function(x) ifelse(x[which(colnames(SNV.data)=="Variant_Classification")]=="indel",
                                                               ifelse(length(grep("frameshift", x[which(colnames(SNV.data)=="Plot_Class")]))>0,
                                                                      ifelse(x[which(colnames(SNV.data)=="Reference_Allele")]=="-",
                                                                             "frameshift_INS",
                                                                             "frameshift_DEL"),
                                                                      ifelse(length(grep("inframe_insertion", x[which(colnames(SNV.data)=="Plot_Class")]))>0,
                                                                             "inframe_insertion",
                                                                             ifelse(length(grep("inframe_deletion", x[which(colnames(SNV.data)==
                                                                                                                              "Plot_Class")]))>0,
                                                                                    "inframe_deletion",
                                                                                    x[which(colnames(SNV.data)=="Plot_Class")]))),
                                                               ifelse(x[which(colnames(SNV.data)=="Plot_Class")]=="initiator_codon",
                                                                      ifelse(x[which(colnames(SNV.data)=="Variant_Classification")]=="Missense_Mutation",
                                                                             "missense",
                                                                             x[which(colnames(SNV.data)=="Plot_Class")]),
                                                                      x[which(colnames(SNV.data)=="Plot_Class")]))
                                            ))
    SNV.data$Plot_Class = gsub("deletion", "DEL", gsub("insertion", "INS", SNV.data$Plot_Class))
  } else {
    SNV.data %<>% mutate(Plot_Class = apply(SNV.data, 1,
                                            function(x) ifelse(x[which(colnames(SNV.data)=="Variant_Classification")]=="RNA" &&
                                                                 length(grep("splice", x[which(colnames(SNV.data)=="Effect_Ontology")]))>0 &&
                                                                 length(grep("SPLICE", x[which(colnames(SNV.data)=="Effect_Class")]))>0,
                                                               "Splice_Site",
                                                               x[which(colnames(SNV.data)=="Variant_Classification")])))
    SNV.data %<>% mutate(Plot_Class = apply(SNV.data, 1,
                                            function(x) ifelse(x[which(colnames(SNV.data)=="Plot_Class")]=="indel",
                                                               ifelse(length(grep("frameshift", x[which(colnames(SNV.data)=="Effect_Ontology")]))>0,
                                                                      ifelse(x[which(colnames(SNV.data)=="Reference_Allele")]=="-",
                                                                             "Frame_Shift_INS",
                                                                             "Frame_Shift_DEL"),
                                                                      ifelse(length(grep("inframe_insertion", x[which(colnames(SNV.data)=="Effect_Ontology")]))>0,
                                                                             "In_Frame_Ins",
                                                                             ifelse(length(grep("inframe_deletion", x[which(colnames(SNV.data)==
                                                                                                                              "Effect_Ontology")]))>0,
                                                                                    "In_Frame_Del",
                                                                                    x[which(colnames(SNV.data)=="Plot_Class")]))),
                                                               x[which(colnames(SNV.data)=="Plot_Class")])))
  }
  SNV.data %<>% filter(grepl(paste(SNV_mutations, collapse="|"), Plot_Class))
  return(list(data = SNV.data, stats = SNV.stats))
}

# ***************
# Function to populate df with cohort SNV 
# ***************
populate.sv = function(study.data,
                       temp.dir=temp_dir
                       ){
  # ***************
  # Set up dataframe
  # ***************
  SV.data = NULL
  sv.table=NULL
  # ***************
  # Populate dataframe
  # ***************
  for (i in 1:nrow(study.data)){
    # ***************
    # read in SV data
    # ***************
    sample.file = ifelse(!is.na(study.data$testCaptureKit[i])&!study.data$testCaptureKit[i]=="", #nchar(study.data$testCaptureKit[i])>0, 
                         paste(study.data$analysisPath[i], SV_WES_file, sep="/"),
                         paste(study.data$qsvAnalysisPath[i], SV_WGS_file, sep="/"))
    if (system(paste("ls -A ", sample.file, " 2>/dev/null | wc -l", sep=""),intern=TRUE)==1){
      # message(paste("Structural variants found for", study.data$UniqueLabel[i]))
      SV.filename = system(paste("ls -A ", sample.file, " 2>/dev/null ", sep=""),intern=TRUE)
      del.header = as.integer(system(paste("cat ", SV.filename, " | grep -c '#' ", sep=""), intern=TRUE))
      sample.data = fread(SV.filename, skip=del.header,
                          na.strings="NA", strip.white=TRUE,
                          colClasses=c("character"), stringsAsFactors=FALSE)
      # add a line here to accomodate the BGI/Illumina WGS sequencing differences in the SCLC WGS analysis
      sample.data %<>% dplyr::filter(category=="1"|category=="2")
      sample.table = data.frame(UniqueLabel = c(study.data$UniqueLabel[i]),
                                FilePath = ifelse(!is.na(study.data$testCaptureKit[i]),
                                                  study.data$analysisPath[i],
                                                  study.data$qsvAnalysisPath[i]),
                                FileExists = system(paste("ls -A ", sample.file, " 2>/dev/null | wc -l", sep=""),intern=TRUE)==1,
                                SV_num = nrow(sample.data))
      if (is.null(sv.table)){
        sv.table = sample.table
      } else {
        sv.table = rbind(sv.table, sample.table)
      }
      
    } else {
      message(paste("Structural variants", ifelse(length(grep("WES", study.data$UniqueLabel[i]))>0, "", sample.file), 
                    "NOT found for", study.data$UniqueLabel[i]))
      sample.data = NULL
    }
    # ***************
    # create SV_data
    # ***************
    if (!is.null(sample.data) && nrow(sample.data) != 0){
      message(paste("SV data has been read in for ", study.data$UniqueLabel[i], sep=""))
      sample.data %<>% mutate(DonorLabel = study.data$donorLabel[i],
                              UniqSampleLabel = study.data$UniqueLabel[i]) %>% select(DonorLabel,
                                                                                      UniqSampleLabel,
                                                                                      chr_from,
                                                                                      chr_from_bkpt,
                                                                                      chr_from_strand,
                                                                                      chr_from_range,
                                                                                      gene_name_from,
                                                                                      bkpt_from_context,
                                                                                      chr_to,
                                                                                      chr_to_bkpt,
                                                                                      chr_to_strand,
                                                                                      chr_to_range,
                                                                                      gene_name_to,
                                                                                      bkpt_to_context,
                                                                                      annotation,
                                                                                      type,
                                                                                      likely_consequence)
      SV.data = rbind(SV.data, sample.data)
    }
  }
  return(list(data = SV.data, stats=sv.table))
}

# ***************
# Function to populate df with cohort SNV 
# ***************
populate.cnv = function(study.data,
                        samp.type="WGS",
                        temp.dir=temp_dir,
                        gene.list = Gene_list){
  # ***************
  # Set up dataframe
  # ***************
  CNV.data = NULL
  # ***************
  # Populate dataframe
  # ***************
  if (samp.type=="WGS"){
    cnv.type = "ascat"
  } else if (samp.type=="WES") {
    cnv.type = "sequenza"
  }
  for (i in 1:nrow(study.data)){
    # ***************
    # read in CN data
    # ***************
    if (!is.na(study.data$ascatAnalysisPath[i]) && (is.na(study.data$testCaptureKit[i])||study.data$testCaptureKit[i]=="")){
      message(paste("Ascat analysis was run for", study.data$UniqueLabel[i]))
      #
      # nominate where ascat results for WGS are stored
      #
      wgs.ascat.path = ""
      sample.file = list.files(path=paste0(wgs.ascat.path, "cnv_data/", study.data$donorLabel[i]), pattern=CNV_file, full.names=TRUE)
      ascat.file = list.files(path=paste0(wgs.ascat.path, "cnv_data/", study.data$donorLabel[i]), pattern=".copynumber.ascat.tsv", full.names=TRUE)
      # sample.file = paste(study.data$ascatAnalysisPath[i], "/",
      #                    unlist(strsplit(study.data$ascatAnalysis[i], ":"))[2], CNV_file, sep="")
    } else if (!is.na(study.data$sequenzaAnalysisPath[i]) && length(grep("exo|nextera", study.data$testCaptureKit[i], ignore.case=T))>0){
      message(paste("Sequenza analysis was run for", study.data$UniqueLabel[i]))
      sample.file = list.files(path=study.data$sequenzaAnalysisPath[i], pattern=CNV_file, full.names=TRUE)
      # sample.file = paste(study.data$sequenzaAnalysisPath[i], "/",
      #                     unlist(strsplit(study.data$sequenzaAnalysis[i], ":"))[2], CNV_file, sep="")
    } else {
      message(paste("No ascat or sequenza CNV analysis found for", study.data$UniqueLabel[i]))
      sample.file = ""
      sample.data = NULL
      cnv.type = "none"
    }
    if (!sample.file==""){
      if (file.exists(sample.file)) {
        # message(paste("Read in", cnv_type, "CNV file found for", study.data$UniqueLabel[i]))
        sample.data = as.data.frame(fread(sample.file, header=TRUE, na.strings = "NA",
                                          strip.white=TRUE, colClasses=c("character"), stringsAsFactors=FALSE))
        sample.data %<>% filter(!is.na(Copynumber))
        if (cnv.type=="ascat"){
          sample.data %<>% mutate(Copy_Number_Variant = ASCAT_call)
          if (file.exists(ascat.file)){
            ascat.data = as.data.frame(fread(ascat.file, header=TRUE, na.strings = "NA",
                                             strip.white=TRUE, colClasses=c("character"), stringsAsFactors=FALSE))
            ascat.data %<>% filter(tumour.minor.cn==0)
          }
        } else if (cnv.type=="sequenza"){
          sample.data %<>% mutate(Copy_Number_Variant = Sequenza_call)
        }
        sample.file = ""
      } else {
        message(paste(cnv.type, "Copy Number file", sample.file, "NOT found for", study.data$UniqueLabel[i]))
        sample.file = ""
        sample.data = NULL
      }
    } else {
      message(paste("No", cnv.type, "CNV analysis found for", study.data$UniqueLabel[i]))
      sample.file = ""
      sample.data = NULL
    }
    # ***************
    # create CNV_data
    # ***************
    # limit CNV sample data to the genes of interest
    if (!is.null(sample.data) && !is.null(nrow(sample.data)) && nrow(sample.data) > 0){
      # message(paste("CNV data has been found for ", study.data$UniqueLabel[i], sep=""))
      sample.data %<>% mutate(DonorLabel = study.data$donorLabel[i],
                              UniqSampleLabel = study.data$UniqueLabel[i],
                              Hugo_Symbol = gene_symbol) %>% select(DonorLabel,
                                                                    UniqSampleLabel,
                                                                    Hugo_Symbol,
                                                                    Copy_Number_Variant,
                                                                    Copynumber,
                                                                    chromosome,
                                                                    gene_start,
                                                                    gene_end,
                                                                    SegmentDescription,
                                                                    PercentGeneCovered)
      if (cnv.type=="ascat"){
        for (j in 1:nrow(sample.data)){
          ascat.chr = dplyr::filter(ascat.data, chromosome==sample.data$chromosome[j], !event=="HET")
          if (sample.data$SegmentDescription[j]=="span"){
            ascat.chr %<>% dplyr::filter(start<=sample.data$gene_start[j], end>=sample.data$gene_end[j])
          } else if (sample.data$SegmentDescription[j]=="partial") {
            cond1 = dplyr::filter(ascat.chr, start>sample.data$gene_start[j]&start<sample.data$gene_end[j]&end>sample.data$gene_end[j])
            cond2 = dplyr::filter(ascat.chr, start<sample.data$gene_start[j]&end<sample.data$gene_end[j]&end>sample.data$gene_start[j])
            ascat.chr = rbind(cond1, cond2)
            if (nrow(ascat.chr)>1){
              
            }
          } else {
            ascat.chr %<>% dplyr::filter(start>sample.data$gene_start[j]&end<sample.data$gene_end[j])
          }
          sample.data$LOH[j] = ifelse(nrow(ascat.chr)>0, paste(unique(ascat.chr$LOH), collapse=";"), "FALSE")
        }
      }
      CNV.data = rbind(CNV.data, sample.data)
    } else {
      message(paste("No CNV data was found for", study.data$UniqueLabel))
    }
  } 
  return(CNV.data)
}

# ***************
# Function to extract sample cellularity from sequenza and ascat analyses
# ***************
get.cellularity = function(study.data,
                           uniq.label){
  sample.loc = grep(uniq.label, study.data$UniqueLabel)
  if (length(sample.loc)>1){
    if (length(which(!is.na(study.data$ascatAnalysisPath[sample.loc])))>0){
      cell.loc = which(!is.na(study.data$ascatAnalysisPath[sample.loc]))
    } else if (length(which(!is.na(study.data$sequenzaAnalysisPath[sample.loc])))>0){
      cell.loc = which(!is.na(study.data$sequenzaAnalysisPath[sample.loc]))
    } else {
      cell.loc = 1
    }
    sample.loc = sample.loc[cell.loc[1]]
  }
  if (length(sample.loc)>0){
    if (!study.data$ascatAnalysisPath[sample.loc]==""||!study.data$sequenzaAnalysisPath[sample.loc]==""){
      if (!study.data$ascatAnalysisPath[sample.loc]=="") {
        file = paste(study.data$ascatAnalysisPath[sample.loc], "/",
                     gsub("analysis:", "", study.data$ascatAnalysis[sample.loc]),".summary.stats.tsv", sep="")
        purity_var = "ascat"
      } else {
        file = paste(study.data$sequenzaAnalysisPath[sample.loc], "/",
                     gsub("collectedsample:", "", study.data$testSampleIRI[sample.loc]),
                     "_confints_CP.txt", sep="")
        purity_var = "sequenza"
      }
      if (file.exists(file)){
        purity.data = read.delim(file,
                                 header=TRUE, na.strings="NA", strip.white=TRUE,
                                 colClasses=c("character"), stringsAsFactors=FALSE)
        cellularity = ifelse(purity_var == "ascat",
                             as.numeric(purity.data$purity)*100,
                             ifelse(purity_var == "sequenza",
                                    as.numeric(purity.data$cellularity[2])*100, "NA"))
      } else {
        cellularity = "NA"
      }
    } else {
      cellularity = "NA"
    }
    cellularity = ifelse(is.na(cellularity), "NA", cellularity)
  } else {
    cellularity = "NA"
  } 
  return(cellularity)
}

# ***************
# Function to create cellularity status
# ***************
get.cellularity.status = function(study.data,
                                  uniq.label){
  cellularity = study.data$Cellularity[which(study.data$donorLabel==uniq.label)]
  if (!is.na(cellularity)){
    cellularity.status = ifelse(cellularity=="NA",
                                "NA",
                                ifelse(cellularity>=90, ">=90%",
                                       ifelse(cellularity>=75, "75% to <90%",
                                              ifelse(cellularity>=50, "50% to <75%",
                                                     ifelse(cellularity<50, "<50%", "NA")))))
  } else {
    cellularity.status = "NA"
  }
  cellularity.status = paste("Cell", cellularity.status, sep="_")
  return(cellularity.status)
}

# ***************
# Function to extract sample ploidy from ascat and sequenza analyses
# ***************
get.ploidy = function(study.data,
                      uniq.label){
  sample.loc = which(study.data$UniqueLabel==uniq.label)
  if (!is.na(study.data$ascatPloidyScore[sample.loc])){
    ploidy = study.data$ascatPloidyScore[sample.loc]
  } else if (!is.na(study.data$sequenzaPloidyScore[sample.loc])){
    ploidy = study.data$sequenzaPloidyScore[sample.loc]  
  } else {
    ploidy = "NA"
  }
  return(ploidy)
}

# ***************
# Function to create sample tumour mutation burden status
# ***************
get.tmb.status = function(study.data,
                          uniq.label){
  tmb = study.data$Adj_TMB[which(study.data$UniqueLabel==uniq.label)]
  if (!is.na(tmb)){
    tmb.status = ifelse(tmb>15, ">15",
                        ifelse(tmb>10, ">10 to <=15",
                               ifelse(tmb>5, ">5 to <=10", 
                                      ifelse(tmb>=0, "<=5", "NA"))))
  } else {
    tmb.status = "NA"
  }
  tmb.status = paste("TMB", tmb.status, sep="_")
  return(tmb.status)
}

# ***************
# Function to create sample microsatellite instability status
# ***************
get.msi.status = function(study.data,
                          uniq.label){
  msi = study.data$MSI_perc[which(study.data$UniqueLabel==uniq.label)]
  if (!is.na(msi)){
    msi.status = ifelse(msi>=3.5, ">3.5",
                        ifelse(msi>=0, "<3.5", "NA"))
  } else {
    msi.status = "NA"
  }
  msi.status = paste("MSI", msi.status, sep="_")
  return(msi.status)
}

# *************
# Function to add the main annotations to the study data
# *************
add.annotations = function(study.data=study_data,
                           top.annots,
                           study.prefix = study_prefix,
                           work.dir = working_dir,
                           case = "UniqueLabel",
                           status = TRUE,
                           status.file = Status_data,
                           subgp.col = subgp_colname,
                           samp.id = samp_id){
  capture.stats = fread(paste(gsub("GeneLandscape_Summary", "Extra_data", work.dir), "Fielding_-_Lung_Illumina_357_WGS_WES_capture_Stats.tsv", sep=""))
  if ("TMB" %in% top.annots){
    study.data$Adj_TMB = sapply(study.data$UniqueLabel, function(x) ifelse(length(grep("WGS", x))>0, 
                                                                             ifelse("WGS_UniqLabel" %in% colnames(status.file), status.file$WGS_Nonsyn_TMB[which(status.file$WGS_UniqLabel==x)], status.file$WGS_Nonsyn_TMB[which(status.file$UniqueLabel==x)]),
                                                                             ifelse("WES_UniqLabel" %in% colnames(status.file), status.file$WES_Nonsyn_TMB[which(status.file$WES_UniqLabel==x)], status.file$WES_Nonsyn_TMB[which(status.file$UniqueLabel==x)])))
    study.data %<>% dplyr::mutate(TMB = sapply(as.character(study.data[,case]), function(x) get.tmb.status(study.data, x), USE.NAMES=FALSE))
  }
  if ("MSI" %in% top.annots){
    study.data %<>% mutate(MSI_perc = sapply(study.data[,case], function(x) ifelse(x%in%capture.stats$UniqueLabel,
                                                                                   capture.stats$MSI[which(capture.stats$UniqueLabel==x)],
                                                                                   1*-1), USE.NAMES=FALSE))
    study.data %<>% mutate(MSI = sapply(as.character(study.data[,case]), function(x) get.msi.status(study.data, x), USE.NAMES=FALSE))
  }
  if ("Meth_Cluster" %in% top.annots){
    study.data %<>% mutate(Meth_Cluster=sapply(study.data$UniqueLabel, function(x) status.file$Meth_Cluster[grep(substr(x,1,10), status.file[,samp.id])]))
  }
  if ("WGD_Status" %in% top.annots){
    study.data %<>% mutate(WGD_Status=sapply(study.data$UniqueLabel, function(x) status.file$CN_status[grep(substr(x,1,10), status.file[,samp.id])]))
  }
  if ("Amp_Arch" %in% top.annots){
    study.data %<>% mutate(Amp_Arch=sapply(study.data$UniqueLabel, function(x) status.file$Amp_Arch[grep(substr(x,1,10), status.file[,samp.id])]))
  }
  if ("Treatment" %in% top.annots){
    study.data %<>% mutate(Treatment=sapply(study.data$UniqueLabel, function(x) status.file$Treatment[grep(substr(x,1,10), status.file[,samp.id])]))
  }
  if ("Survival" %in% top.annots){
    study.data %<>% mutate(Survival=sapply(study.data$UniqueLabel, function(x) status.file$Survival[grep(substr(x,1,10), status.file[,samp.id])]))
  }
  if ("Cellularity" %in% top.annots){
    study.data %<>% mutate(Cellularity=sapply(study.data$UniqueLabel, function(x) get.cellularity.status(status.file, substr(x,1,10))))
    # study.data %<>% mutate(Cellularity = sapply(as.character(study.data[,case]), function(x) get.cellularity.status(study.data, x), USE.NAMES=FALSE))
    study.data %<>% arrange(factor(Cellularity, levels=c("Cell_<50%", "Cell_50% to <75%", "Cell_75% to <90%", "Cell_>=90%")), UniqueLabel) 
  }
  if (status){
    # Genomic_reclassification_telomere.csv
    Donor_status = data.frame(Donor_ID = study.data$UniqueLabel,
                              Donor_class = sapply(study.data$UniqueLabel, function(x) ifelse(is.na(unique(status.file[grep(substr(x,1,10), status.file[,samp.id]), subgp.col])),
                                                                                              "other",
                                                                                              unique(status.file[grep(substr(x,1,10), status.file[,samp.id]), subgp.col]))))
   
    # set levels 
    Donor_status$Donor_class = factor(Donor_status$Donor_class, levels=c("adeno", "NSCLC", "squamous", "sclc", "LC_neuro", "other_met", "other"))
    Donor_status$Donor_class = droplevels(Donor_status$Donor_class)
    study.data = cbind(study.data, Donor_status[match(study.data$UniqueLabel, Donor_status$Donor_ID),])
    study.data %<>% arrange(Donor_class, UniqueLabel)
  } 
  return(list(study=study.data))
}

# ***************
# Function to extract and contract the effect of edited mutation types 
# ***************
# Is designed to report the effect of highest consequence between Nonsense, Frameshift or Missense.
get.effect = function(full.effect, 
                      effect,
                      type,
                      vep.anno=VEP_anno){
  if (vep.anno) {
    new.effect = switch(effect,
                        stop_gained = ifelse(length(grep("frameshift", full.effect))>0, 
                                             gsub("frameshift", "stop_gained", full.effect),
                                             ifelse(length(grep("missense", full.effect))>0, 
                                                    gsub("missense", "stop_gained", full.effect), 
                                                    paste(effect, full.effect, sep=";"))),
                        stop_lost = ifelse(length(grep("frameshift", full.effect))>0, 
                                           gsub("frameshift", "stop_lost", full.effect),
                                           ifelse(length(grep("missense", full.effect))>0, 
                                                  gsub("missense", "stop_lost", full.effect), 
                                                  paste(effect, full.effect, sep=";"))),
                        frameshift= ifelse(length(grep("stop_gained", full.effect))>0, 
                                           full.effect,
                                           ifelse(length(grep("splice", full.effect))>0,
                                                  full.effect,
                                                  ifelse(length(grep("missense", full.effect))>0, 
                                                         gsub("missense", "frameshift", full.effect),
                                                         paste(effect, full.effect, sep=";")))),
                        missense = ifelse(length(grep("stop_gained|splice_acceptor|splice_donor|frameshift", full.effect))>0, 
                                          full.effect, 
                                          paste(effect, full.effect, sep=";")),
                        new = full.effect)
  } else {
    new.effect = switch(effect,
                        Nonsense_Mutation = ifelse(length(grep("Frame_Shift", full.effect))>0, 
                                                   gsub("Frame_Shift", "Nonsense_Mutation", full.effect),
                                                   ifelse(length(grep("Missense_Mutation", full.effect))>0, 
                                                          gsub("Missense_Mutation", "Nonsense_Mutation", full.effect), 
                                                          paste(effect, full.effect, sep=";"))),
                        Frame_Shift = ifelse(length(grep("Nonsense_Mutation", full.effect))>0, 
                                             full.effect,
                                             ifelse(length(grep("Missense_Mutation", full.effect))>0, 
                                                    gsub("Missense_Mutation", "Frame_Shift", full.effect),
                                                    paste(effect, full.effect, sep=";"))),
                        Missense_Mutation = ifelse(length(grep("Nonsense_Mutation|Frame_Shift", full.effect))>0, 
                                                   full.effect, 
                                                   paste(effect, full.effect, sep=";")),
                        new = full.effect)
  }
  other.variants = c("transcription_ablation",  "splice_acceptor", "splice_donor", "splice_region", "splice_site",
                     "start_gained", "start_lost", "transcription_amplification", "inframe",
                     "protein_altering", "stop_retained", "start_retained", "incomplete_terminal_codon", 
                     "regulatory_region_ablation", 
                     "In_Frame", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")
  if (length(grep(effect, other.variants))>0){
    new.effect = paste(effect, paste(full.effect, collapse=";"), sep=";")
  }
  new.effect = unlist(strsplit(new.effect, ";"))
  return(new.effect)
}


# ***************
# Function to extract and contract the effect of edited mutation types 
# ***************
# Is designed to report the effect of highest consequence between Nonsense, Frameshift or Missense.
get.most.del = function(gene.effect,
                        comp.effect = NULL,
                        top.genes = top_genes,
                        add.ctdna = add_ctdna){
  gp1.var = c("stop_gained", "stop_lost", "start_gained", "start_lost", "frameshift")
  gp2.var = c("missense", "inframe")
  other.var = c("transcription_ablation",  "splice_acceptor", "splice_donor", "splice_region", "splice_site",
                "transcription_amplification", "protein_altering", "stop_retained", "start_retained", "incomplete_terminal_codon", 
                "regulatory_region_ablation")
  if (!is.null(comp.effect)){
    if (length(comp.effect)>1){
      if (length(which(comp.effect%in%gp1.var))>0){
        most.del = comp.effect[which(comp.effect%in%gp1.var)]
      } else if (length(which(comp.effect%in%gp2.var))>0){
        most.del = comp.effect[which(comp.effect%in%gp2.var)]
      } else if (length(which(comp.effect%in%other.var))>0){
        most.del = comp.effect[which(comp.effect%in%other.var)]
      }
    } else {
      most.del = comp.effect
    }
    gp1.var = unique(c(most.del, gp1.var))
    gp2.var = gp2.var[which(!gp2.var%in%most.del)]
    other.var = other.var[which(!other.var%in%most.del)]
  }
  gp1.effect = gene.effect[which(gene.effect%in%gp1.var)]
  gp2.effect = gene.effect[which(gene.effect%in%gp2.var)]
  other.effect = gene.effect[which(gene.effect%in%other.var)]
  gp1 = ifelse(length(gp1.effect)>0, TRUE, FALSE)
  gp2 = ifelse(length(gp2.effect)>0, TRUE, FALSE)
  other = ifelse(length(other.effect)>0, TRUE, FALSE)
  new.effect = ifelse((gp1&gp2)|(gp1&other),
                      ifelse(!is.null(top.genes)|add.ctdna, paste0(gp1.effect[order(gp1.effect, decreasing = T)][1], ";+"), "multi-hit"),
                      ifelse(gp1,
                             ifelse(length(gp1.effect)>length(unique(gp1.effect)),
                                    ifelse(!is.null(top.genes)|add.ctdna, paste0(gp1.effect[order(gp1.effect, decreasing = T)][1], ";+"), "multi-hit"),
                                    gp1.effect[order(gp1.effect, decreasing = T)][1]),
                             ifelse(gp2&other,
                                    ifelse(!is.null(top.genes)|add.ctdna, paste0(gp2.effect[order(gp2.effect, decreasing = T)][1], ";+"), "multi-hit"),
                                    ifelse(gp2,
                                           ifelse(length(gp2.effect)>length(unique(gp2.effect)),
                                                  ifelse(!is.null(top.genes)|add.ctdna, paste0(gp2.effect[order(gp2.effect, decreasing = T)][1], ";+"), "multi-hit"),
                                                  gp2.effect[order(gp2.effect, decreasing = T)][1]),
                                           ifelse(length(other.effect)>length(unique(other.effect)),
                                                  ifelse(!is.null(top.genes)|add.ctdna, paste0(other.effect[order(other.effect, decreasing = T)][1], ";+"), "multi-hit"),
                                                  other.effect[order(other.effect, decreasing = T)][1])))))
  return(new.effect)
}
# ***************
# Function to determine copy number variant to report in summary.
# ***************
get.copy = function(effect, 
                    copynum, 
                    loh,
                    type, 
                    flag.loss, 
                    flag.gain,
                    study.data,
                    uniq.label){
  ploidy = get.ploidy(study_data,
                      as.character(uniq.label))
  # editted to report either HighGain or Gain depending on ploidy and flag
  # editted to report TotalLoss if CN0 or if ploidy >=2.7 CN-2.7
  switch(effect,
         Del = ifelse(type=="all", 
                      "LOH", # "Loss", 
                      ifelse(copynum==0, "TotalLoss", "LOH")), #ifelse(copynum==0, "TotalLoss", "Loss")), #type=="all", "Loss", ifelse(
         cnLOH = ifelse(flag.loss, "LOH","cn2LOH"), 
         Gain = ifelse(flag.gain, 
                       ifelse(copynum>=9, 
                              "HighGain", 
                              ifelse(copynum>=5, 
                                     ifelse(ploidy<=2.7, 
                                            "HighGain", 
                                            ifelse(loh==TRUE || loh=="TRUE", "LOH", "DontReport")),
                                     ifelse(loh==TRUE || loh=="TRUE", "LOH", "DontReport"))), 
                       ifelse(copynum>=9, 
                              "HighGain", 
                              ifelse(copynum>=5, 
                                     ifelse(ploidy<=2.7, 
                                            "HighGain", 
                                            ifelse(loh==TRUE || loh=="TRUE", "LOH", "Gain")),
                                     ifelse(loh==TRUE || loh=="TRUE", "LOH", "DontReport")))),
         Het = "DontReport",
         HighGain = ifelse(flag.gain, 
                           ifelse(copynum>=9, 
                                  "HighGain", 
                                  ifelse(copynum>=5, 
                                         ifelse(ploidy<=2.7, 
                                                "HighGain", 
                                                ifelse(loh==TRUE || loh=="TRUE", "LOH", "DontReport")),
                                         ifelse(loh==TRUE || loh=="TRUE", "LOH", "DontReport"))), 
                           ifelse(copynum>=9, 
                                  "HighGain", 
                                  ifelse(copynum>=5, 
                                         ifelse(ploidy<=2.7, 
                                                "HighGain", 
                                                ifelse(loh==TRUE || loh=="TRUE", "LOH", "Gain")),
                                         ifelse(loh==TRUE || loh=="TRUE", "LOH", "DontReport"))))
  )
}

# ***************
# Function to determine the new summarised variant to report in summary.
# ***************
get.new.effect = function(SNV.df=NULL, 
                          CNV.df=NULL, 
                          SV.df=NULL, 
                          comp.data = NULL, # for ctDNA comparison
                          type,
                          study.data,
                          uniq.label){
  new.effect = NULL
  gp1.var = c("stop_gained", "stop_lost", "start_gained", "start_lost", "frameshift_INS", "frameshift_DEL", "frameshift")
  gp2.var = c("missense", "inframe", "inframe_INS", "inframe_DEL", "inframe_insertion", "inframe_deletion")
  other.var = c("transcription_ablation",  "splice_acceptor", "splice_donor", "splice_region", "splice_site",
                "transcription_amplification", "protein_altering", "stop_retained", "start_retained", "incomplete_terminal_codon", 
                "regulatory_region_ablation")
  if (type %in% c("SNV", "all")){
    if (!is.null(SNV.df) && nrow(SNV.df)>0){
      gene.effect = NULL
      comp.effect = NULL
      for (i in 1:nrow(SNV.df)){
        effect = unlist(strsplit(as.character(SNV.df$Plot_Class[i]), ";"))
        if (length(effect)>1){
          effect = ifelse(length(which(effect%in%gp1.var))>0,
                          effect[which(effect%in%gp1.var)],
                          ifelse(length(which(effect%in%gp2.var))>0,
                                 effect[which(effect%in%gp2.var)],
                                 ifelse(length(which(effect%in%other.var))>0,
                                        effect[which(effect%in%other.var)],
                                        paste(effect, collapse=";"))))
        }
        # report Frame shift insertions or deletions as Frame shifts
        effect = ifelse(length(grep("Frame_Shift|frameshift", effect))>0,
                        ifelse(length(grep("Frame_Shift", effect))>0,
                               "Frame_Shift",
                               "frameshift"),
                        # report inframe insertions or deletions as inframe
                        ifelse(length(grep("In_Frame|inframe", effect))>0,
                               ifelse(length(grep("In_Frame", effect))>0,
                                      "In_Frame",
                                      "inframe"),
                               # report splice site acceptor, donor, or region as splice_stie
                               ifelse(length(grep("splice", effect))>0,
                                      "splice_site",
                                      effect)))
        if (!is.null(comp.data)){
          wgs.effect = dplyr::filter(comp.data, VarID%in%SNV.df$VarID[i])$Plot_Class
          if (length(wgs.effect)>0){
            wgs.effect = ifelse(length(grep("Frame_Shift|frameshift", wgs.effect))>0,
                                ifelse(length(grep("Frame_Shift", wgs.effect))>0,
                                       "Frame_Shift",
                                       "frameshift"),
                                # report inframe insertions or deletions as inframe
                                ifelse(length(grep("In_Frame|inframe", wgs.effect))>0,
                                       ifelse(length(grep("In_Frame", wgs.effect))>0,
                                              "In_Frame",
                                              "inframe"),
                                       # report splice site acceptor, donor, or region as splice_stie
                                       ifelse(length(grep("splice", wgs.effect))>0,
                                              "splice_site",
                                              wgs.effect)))
            comp.effect = c(comp.effect, wgs.effect)
            if (!wgs.effect==effect){
              if (length(grep(wgs.effect, gp1.var))>0){
                effect = wgs.effect
              } else if (length(grep(wgs.effect, gp2.var))>0 & !length(grep(effect, gp1.var))>0){
                effect = wgs.effect
              }
            }
          }
        }
        if (effect=="RNA"){ 
          gene.effect = gene.effect 
        } else {
          gene.effect = c(gene.effect, effect)
        }
      }
      new.effect = get.most.del(gene.effect = gene.effect, comp.effect = comp.effect)
    }
  }
  if (type %in% c("all")){
    if (!is.null(CNV.df) && nrow(CNV.df)>0){
      for (k in 1:nrow(CNV.df)){
        effect = as.character(CNV.df$Copy_Number_Variant[k])
        effect.before = effect
        copynum = as.integer(CNV.df$Copynumber[k])
        loh = CNV.df$LOH[k]
        effect = get.copy(effect=effect, 
                          copynum=as.integer(copynum), 
                          loh = loh,
                          type=type, 
                          flag.loss=report_just_loss, 
                          flag.gain=report_just_HG,
                          study.data,
                          uniq.label)
        new.effect = unlist(ifelse(effect=="DontReport", 
                                   list(new.effect),
                                   ifelse(is.null(new.effect) || length(new.effect) == 0, 
                                          effect,
                                          ifelse(length(grep(effect, new.effect)) == 0, 
                                                 paste(new.effect, effect, sep=";"), 
                                                 new.effect))))
      }
    }
  }
  if (type %in% c("all")){
    if (!is.null(SV.df) && nrow(SV.df)>0){
      for (j in 1:nrow(SV.df)){
        if (type == "all"){
          if (length(grep("Loss_of_Function|Fusion", SV.df$likely_consequence[j], ignore.case=TRUE))>0){
            new.effect = ifelse(is.null(new.effect) || length(new.effect) == 0, 
                                "SV", 
                                ifelse(length(grep("SV", new.effect)) == 0,
                                       paste(new.effect, "SV", sep=";"), new.effect))
          }
        } else {
          effect = SV.df$Reported_consequence[j]
          if (is.null(new.effect) || length(new.effect) == 0){
            new.effect = effect
          } else if (length(grep(effect, new.effect)) == 0) {
            new.effect = c(new.effect, effect)
          } else {
            new.effect = new.effect
          }
        }
      }
    }
  }
  new.effect = ifelse(is.null(new.effect), NA, new.effect)
  return(new.effect)
}

# ***************
# Function to extract and sort the data to be plotted
# ***************
prepare.data = function(df, 
                        type, 
                        mut.types, 
                        variant.class,
                        gene.list = Gene_list){
  if (!is.null(df) && nrow(df) != 0){
    #extract the data to be plotted
    if (type == "SV"){
      df$gene_name_from = unlist(lapply(strsplit(df$gene_name_from, "|", fixed=TRUE), function(x) ifelse(length(which(x%in%gene.list$Hugo_Symbol))>0, x[which(x%in%gene.list$Hugo_Symbol)], x)))
      df$gene_name_from = unlist(lapply(strsplit(df$gene_name_from, ";", fixed=TRUE), function(x) ifelse(length(which(x%in%gene.list$Hugo_Symbol))>0, x[which(x%in%gene.list$Hugo_Symbol)], x)))
      df$gene_name_from = unlist(lapply(strsplit(df$gene_name_from, "/", fixed=TRUE), function(x) ifelse(length(which(x%in%gene.list$Hugo_Symbol))>0, x[which(x%in%gene.list$Hugo_Symbol)], x)))
      df$gene_name_from = unlist(lapply(strsplit(df$gene_name_from, "-", fixed=TRUE), function(x) ifelse(length(which(x%in%gene.list$Hugo_Symbol))>0, x[which(x%in%gene.list$Hugo_Symbol)], x)))
      df$gene_name_to = unlist(lapply(strsplit(df$gene_name_to, "|", fixed=TRUE), function(x) ifelse(length(which(x%in%gene.list$Hugo_Symbol))>0, x[which(x%in%gene.list$Hugo_Symbol)], x)))
      df$gene_name_to = unlist(lapply(strsplit(df$gene_name_to, ";", fixed=TRUE), function(x) ifelse(length(which(x%in%gene.list$Hugo_Symbol))>0, x[which(x%in%gene.list$Hugo_Symbol)], x)))
      df$gene_name_to = unlist(lapply(strsplit(df$gene_name_to, "/", fixed=TRUE), function(x) ifelse(length(which(x%in%gene.list$Hugo_Symbol))>0, x[which(x%in%gene.list$Hugo_Symbol)], x)))
      df$gene_name_to = unlist(lapply(strsplit(df$gene_name_to, "-", fixed=TRUE), function(x) ifelse(length(which(x%in%gene.list$Hugo_Symbol))>0, x[which(x%in%gene.list$Hugo_Symbol)], x)))
      df = df[union(which(df$gene_name_from%in%as.character(gene.list$Hugo_Symbol)),
                    which(df$gene_name_to%in%as.character(gene.list$Hugo_Symbol))),]
    } else {
      df = df[which(df$Hugo_Symbol %in% as.character(gene.list$Hugo_Symbol)),]
    }
    if (!nrow(df)>0){
      df = NULL
    }
  } else {
    df = NULL
  }
  return(df)
}

# ***************
# Function to create the mutation landscape matrices to report the genes of interest
# ***************
create.mut.landscape = function(type, 
                                df.list,
                                donor.list=Donor_list,
                                gene.list = Gene_list,
                                study.data,
                                top.annos = c("Ploidy", "Cellularity"),
                                top.genes=NULL,
                                add.ctdna=FALSE){
  if (!type%in%c("SNV", "all", "TOP")){
    stop("type paramater must be one of SNV for single nucleotides variants, all for SNV, CNV & SV variant types or TOP to display modified annotations of all variant types with extra information of entered list of genes")
  }
  if (type%in%c("SNV")){
    if (length(df.list)==1 & !"Plot_Class"%in%colnames(df.list[[1]])){
      stop("You must supply Plot_Class of variants for types SNV")
    }
  }
  if (type %in% c("all", "TOP")){
    if (length(df.list)>=3){
      for (i in 1:length(df.list)){
        names(df.list)[i] = ifelse("Plot_Class" %in% colnames(df.list[[i]]), 
                                   "SNV", 
                                   ifelse("Reported_consequence" %in% colnames(df.list[[i]]), 
                                          "SV", 
                                          ifelse("Copy_Number_Variant" %in% colnames(df.list[[i]]), 
                                                 "CNV", 
                                                 ifelse("FC" %in% colnames(df.list[[i]]), "CNV", ""))))
      } 
      if (!"SNV"%in%names(df.list) & !"CNV"%in%names(df.list) & !"SV"%in%names(df.list)){
        stop("You must supply a list of three data frames for types all and TOP with Plot_Class of variants for types SNV, Copy_Number_Variant of variants for type CNV and Reported_consequence of variants for type SV")
      }
    } else {
      stop("You must supply a list of three data frames for types all and TOP with Plot_Class of variants for types SNV, Copy_Number_Variant of variants for type CNV and Reported_consequence of variants for type SV")
    }
  } 
  # define mutations matrix
  df.landscape = matrix(nrow=nrow(gene.list), ncol=length(donor.list))
  # populate the column names and the row names 
  rownames(df.landscape) = as.character(gene.list$Hugo_Symbol)
  colnames(df.landscape) = as.character(donor.list)
  # populate the matrices with the effects
  df.effects = NULL
  for (donor in as.character(donor.list)){
    # message(paste("Donor:", donor))
    sample.genes = NULL
    sample.df = list()
    for (data.set in names(df.list)[which(nchar(names(df.list))>1)]){
      # message(paste("Data set:", data.set))
      if (!is.null(df.list[[data.set]])){
        sample.df[[data.set]] = filter(df.list[[data.set]], UniqSampleLabel==donor)
        if (data.set=="SV"){
          sample.genes = unique(c(sample.genes,
                                  filter(sample.df[[data.set]], gene_name_from%in% gene.list$Hugo_Symbol)$gene_name_from,
                                  filter(sample.df[[data.set]], gene_name_to%in% gene.list$Hugo_Symbol)$gene_name_to))
        } else {
          sample.genes = unique(c(sample.genes, sample.df[[data.set]]$Hugo_Symbol))
        }
      }
    }
    sample.gene.list = sample.genes
    #sample.gene.list = union(unlist(sample.genes[1]), union(unlist(sample.genes[2]), unlist(sample.genes[3])))
    if (length(sample.gene.list)>0){
      for (gene in sample.gene.list){
        # message(paste("Gene:", gene))
        if ("SNV" %in% names(sample.df)){
          sample.SNV.df = sample.df$SNV
          gene.SNV.df = unique(filter(sample.SNV.df, Hugo_Symbol==gene))
          if (!nrow(gene.SNV.df)>0){
            gene.SNV.df = NULL
          }
        } else {
          gene.SNV.df = NULL
        }
        if ("CNV" %in% names(sample.df)){
          sample.CNV.df = sample.df$CNV
          gene.CNV.df = unique(filter(sample.CNV.df, Hugo_Symbol==gene))
          if (!nrow(gene.CNV.df)>0){
            gene.CNV.df = NULL
          }
        } else {
          gene.CNV.df = NULL
        }
        if ("SV" %in% names(sample.df)){
          sample.SV.df = sample.df$SV
          gene.SV.df = dplyr::filter(sample.SV.df, gene_name_from==gene|gene_name_to==gene)
          if (!nrow(gene.SV.df)>0){
            gene.SV.df = NULL
          }
        } else {
          gene.SV.df = NULL
        }
        if (!is.null(gene.SNV.df) || !is.null(gene.CNV.df) || !is.null(gene.SV.df)){
          # message(paste("getting new effect for Donor:", donor, "and gene:", gene))
          # summarise variant effect for each gene per sample
          new.effect = get.new.effect(SNV.df=gene.SNV.df, 
                                      CNV.df=gene.CNV.df, 
                                      SV.df=gene.SV.df, 
                                      type=ifelse(type=="TOP", "all", type),
                                      study.data=study.data,
                                      uniq.label=as.character(study.data$UniqueLabel[which(study.data[,case_to_report]==donor)]))
          # populate matrices with their summarised variant effect
          if (!is.na(new.effect)){
            # message(paste("the new effect is:", new.effect))
            df.landscape[gene, donor] = new.effect
            df.effects = c(df.effects, new.effect)
          }
        }
      }
    }
  }
  if (type=="TOP" & !is.null(top.genes)){
    if (add.ctdna) {
      new.samp.info = as.data.frame(fread(file=paste0(study_dir, "Methylation/WGS_SCLC_Analysis/Sample_Sheets/Fielding_-_Lung_WGS_SCLC_with_MethClust_82SampleSheet_v14.csv")))
      new.samp.info %<>% dplyr::filter(ctDNA=="YES")
      donor.list = donor.list[which(gsub("_WGS_IFF", "", gsub("_WGS_BFF", "", gsub("_WGS_BDQ", "", donor.list)))%in%new.samp.info$Sample_Name)]
      ct.dna = as.data.frame(fread(file=paste0(study_dir, "GeneLandscape_Summary/WGS_SCLC_Analysis/RData/Fielding_-_Lung_TSO_ctDNA_DNA-TMB_Trace_397Samp_Germ_Som_flag_unfiltered.tsv")))
      ct.dna %<>% dplyr::mutate(DonorLabel = gsub("-", "_", ct.dna$DonorLabel)) 
      ct.dna %<>% dplyr::filter(DonorLabel%in%gsub("_WGS", "", gsub("_IFF", "", gsub("_BFF", "", gsub("_BDQ", "", donor.list))))) %>% 
        dplyr::mutate(UniqSampleLabel = sapply(DonorLabel, function(x) ifelse(length(grep(x, donor.list))>0, donor.list[grep(x, donor.list)], x))) %>% dplyr::filter(UniqSampleLabel%in%donor.list, `Germline-variant`==FALSE, GeneName%in%top.genes, VAF>0.002) 
      ct.dna %<>% dplyr::rename(Hugo_Symbol = GeneName,
                                Start_Position = Position,
                                Reference_Allele = RefCall,
                                Tumor_Seq_Allele2 = AltCall,
                                Variant_Type = VariantType,
                                Amino_Acid_Change = ProteinChange,
                                CDS_Change = CDSChange) %>% dplyr::mutate(Plot_Class = gsub("_variant", "", gsub("-", "_", ct.dna$Consequence))) %>% dplyr::select(DonorLabel, UniqSampleLabel, Hugo_Symbol, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, VAF, Variant_Type, Amino_Acid_Change, CDS_Change, Plot_Class)
      # hard code one change of annotation as per IGV
      ct.dna$Plot_Class[which(ct.dna$DonorLabel=="D01_21_043"&ct.dna$Hugo_Symbol=="RB1")] = "frameshift_INS"
      ct.dna %<>% dplyr::mutate(VarID=paste(DonorLabel, gsub("chr", "", Chromosome), Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep=";"))
      df.list[["SNV"]] %<>% dplyr::mutate(VarID=paste(DonorLabel, gsub("chr", "", Chromosome), Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep=";"))
      top.genes.list = dplyr::mutate(ct.dna, Hugo_Symbol = sapply(Hugo_Symbol, function(x) paste(x, "ct", sep="_")))
      top.genes.list = rbind(top.genes.list, 
                             dplyr::filter(df.list[["SNV"]], Hugo_Symbol%in%top.genes, UniqSampleLabel%in%ct.dna$UniqSampleLabel) %>%
                               dplyr::select(DonorLabel, UniqSampleLabel, Hugo_Symbol, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, VAF, Variant_Type, Amino_Acid_Change, CDS_Change, Plot_Class, VarID))
      top.genes.list %<>% dplyr::arrange(DonorLabel, Hugo_Symbol)
      top.genes.list %<>% dplyr::mutate(VarID = paste(DonorLabel, gsub("chr", "", Chromosome), Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep=":"))
      write_csv(top.genes.list, file=paste0(results_dir, "Fielding_-_Lung_", paste(top.genes, collapse="_"), "_SNV_ctDNA.csv"))
    }
    for (gene.name in rownames(df.landscape)){
      message(paste("Processing gene", gene.name))
      if (gene.name %in% top.genes){
        if (add.ctdna){
          df.landscape = rbind(df.landscape, sapply(colnames(df.landscape), function(x) get.new.effect(SNV.df=dplyr::filter(ct.dna, Hugo_Symbol==gene.name, UniqSampleLabel==x), type="SNV", comp.data = df.list[["SNV"]], uniq.label = x)))
          rownames(df.landscape)[nrow(df.landscape)] = paste(gene.name, "ct", sep="_")
        } else {
          df.landscape = rbind(df.landscape, sapply(df.landscape[gene.name,], function(x) ifelse(length(grep("LOH",x))>0, "LOH", NA)))
          rownames(df.landscape)[nrow(df.landscape)] = paste(gene.name, "LOH", sep="_")
          cn.changes = c("Loss", "TotalLoss", "HighGain")
          df.landscape = rbind(df.landscape, sapply(df.landscape[gene.name,], function(x) ifelse(length(grep(paste(cn.changes, collapse="|"),x))>0,
                                                                                                 unlist(strsplit(x, ";"))[grep(paste(cn.changes, collapse="|"),unlist(strsplit(x, ";")))], NA)))
          rownames(df.landscape)[nrow(df.landscape)] = paste(gene.name, "CN", sep="_")
        }
      }
      df.landscape[gene.name,] = sapply(df.landscape[gene.name,], function(x) gsub("LOH", "", gsub("LOH;", "", gsub(";LOH", "", x))))
      df.landscape[gene.name,] = sapply(df.landscape[gene.name,], function(x) gsub("Loss", "", gsub(";Loss", "", gsub("Loss;", "", gsub("HighGain", "", gsub(";HighGain", "", gsub("HighGain;", "", gsub("TotalLoss", "", gsub(";TotalLoss", "", gsub("TotalLoss;", "", x))))))))))
      df.landscape[gene.name,] = sapply(df.landscape[gene.name,], function(x) ifelse(x=="", NA, x))
      if (!add.ctdna){
        df.landscape[gene.name,] = sapply(df.landscape[gene.name,], function(x) ifelse(x=="SV",
                                                                                       "SV",
                                                                                       ifelse(length(grep("frameshift|stop_gained|start_gained", x))>0,
                                                                                              ifelse(length(grep("SV", x))>0, "deleterious;SV", "deleterious"),
                                                                                              ifelse(length(grep("missense|inframe|splice_site",x))>0,
                                                                                                     ifelse(length(grep("SV", x))>0,
                                                                                                            "unknown;SV",
                                                                                                            "unknown"),
                                                                                                     x))))
      }
    }
    if (!is.null(top.genes)){
      df.landscape = df.landscape[grep(paste(top.genes, collapse = "|"), rownames(df.landscape)), donor.list]
    }
    df.effects = unique(unlist(strsplit(apply(df.landscape, 1, function(x) paste(unique(x), collapse=":")), ":")))
    df.effects = df.effects[which(!df.effects=="NA")]
    df.effects = df.effects[which(!df.effects=="")]
    df.effects = df.effects[which(!is.na(df.effects))]
  } else if (type=="TOP" & is.null(top.genes)){
    stop ("You must supply a vector of genes to display the extra information and modified annotations for type TOP")
  }
  if (type=="all"){
    # remove LOH from df.landscape
    for (row in rownames(df.landscape)){
      df.landscape[row,] = gsub("LOH", "", gsub("LOH;", "", gsub(";LOH", "", df.landscape[row,])))
      if (!type=="TOP"){
        df.landscape[row,] = sapply(df.landscape[row,], function(x) ifelse(length(grep(";", x, fixed=TRUE))>0, "multi-hit", x))
      }
      df.landscape[row,] = sapply(df.landscape[row,] , function(x) ifelse(x=="", NA, x))
    }
    df.effects = unique(gsub("LOH", "", gsub("LOH;", "", gsub(";LOH", "", df.effects))))
    df.effects = df.effects[which(!df.effects=="")]
    if (!type=="TOP"){
      df.effects[grep(";", df.effects, fixed=TRUE)] = "multi-hit"
      df.effects = unique(df.effects)
    }
  } 
  # Add annotations to the gene landscape
  for (anno in top.annos){
    if (anno=="Cellularity"){
      df.landscape = rbind(Cellularity = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                           ifelse(length(grep("SOC",x))==0, 
                                                                                                  as.character(study.data$Cellularity[grep(x, study.data$UniqueLabel)]),
                                                                                                  "Cell_NA"),
                                                                                           "Cell_NA"), USE.NAMES=FALSE),
                           df.landscape)
      df.effects = c(df.effects, unique(df.landscape["Cellularity",])) 
    }
    if (anno == "Ploidy"){
      df.landscape = rbind(Ploidy = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                      ifelse(is.numeric(get.ploidy(study_data, x)),
                                                                                             ifelse(get.ploidy(study_data, x)<5,
                                                                                                    paste("PL", round(as.numeric(get.ploidy(study_data, x))), sep=""),
                                                                                                    "PL>5"),
                                                                                             paste("PL", get.ploidy(study_data, x), sep="")),
                                                                                      "NA"),
                                           USE.NAMES=FALSE),
                           df.landscape)
      df.effects = c(df.effects, unique(df.landscape["Ploidy",])) 
    } 
    if (anno == "TMB") {
      df.landscape = rbind(TMB = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                   study.data$TMB[grep(x, study.data$UniqueLabel)], 
                                                                                   "TMB_NA"), USE.NAMES=FALSE),
                           df.landscape)
      df.effects = c(df.effects, unique(df.landscape["TMB",])) 
    }
    if (anno == "MSI") {
      df.landscape = rbind(MSI = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                   study.data$MSI[grep(x, study.data$UniqueLabel)], 
                                                                                   "MSI_NA"), USE.NAMES=FALSE),
                           df.landscape)
      df.effects = c(df.effects, unique(df.landscape["MSI",]))
    }
    if (anno=="Status"){
      if (status){
        df.landscape = rbind(Status = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                        as.character(filter(study.data, UniqueLabel==x)$Donor_class),
                                                                                        "Status_NA"), USE.NAMES=FALSE),
                             df.landscape)
        df.effects = c(df.effects, unique(df.landscape["Status",]))
      }
    }
    if (anno=="Meth_Cluster"){
      df.landscape = rbind(Meth_Cluster = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                            as.character(filter(study.data, UniqueLabel==x)$Meth_Cluster),
                                                                                            "Meth_Cluster_NA"), USE.NAMES=FALSE),
                           df.landscape)
      df.effects = c(df.effects, unique(df.landscape["Meth_Cluster",]))
    }
    if (anno=="WGD_Status"){
      df.landscape = rbind(WGD_Status = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                          as.character(filter(study.data, UniqueLabel==x)$WGD_Status),
                                                                                          "WGD_NA"), USE.NAMES=FALSE),
                           df.landscape)
      df.effects = c(df.effects, unique(df.landscape["WGD_Status",]))
    }
    if (anno=="Amp_Arch"){
      df.landscape = rbind(Amp_Arch = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                        ifelse(is.na(filter(study.data, UniqueLabel==x)$Amp_Arch),
                                                                                               "Amp_Arch_NA",
                                                                                               ifelse(!filter(study.data, UniqueLabel==x)$Amp_Arch=="NA",
                                                                                                      as.character(filter(study.data, UniqueLabel==x)$Amp_Arch),
                                                                                                      "Amp_Arch_NA"))), USE.NAMES=FALSE),
                           df.landscape)
      df.effects = c(df.effects, unique(df.landscape["Amp_Arch",]))
    }
    if (anno=="Treatment"){
      df.landscape = rbind(Treatment = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                         as.character(filter(study.data, UniqueLabel==x)$Treatment),
                                                                                         "Treatment_NA"), USE.NAMES=FALSE),
                           df.landscape)
      df.effects = c(df.effects, unique(df.landscape["Treatment",]))
    }
    if (anno=="Survival"){
      df.landscape = rbind(Survival = sapply(colnames(df.landscape), function(x) ifelse(length(grep(x, study.data$UniqueLabel))>0,
                                                                                        as.character(filter(study.data, UniqueLabel==x)$Survival),
                                                                                        "Survival_NA"), USE.NAMES=FALSE),
                           df.landscape)
      df.effects = c(df.effects, unique(df.landscape["Survival",]))
    }
  }
  df.landscape[which(is.na(df.landscape))] = "NA"
  return(list(landscape=df.landscape, effects=df.effects))
}

# ***************
# Function to create the data and parameters for plotting matrix
# ***************
get.plot.data = function (type, 
                          v.effects, 
                          v.df, 
                          v.filename,
                          data, 
                          donor.list=Donor_list,
                          gene.list = Gene_list,
                          study.data = study_data,
                          plot.type = "",
                          gene.order=NULL,
                          results.dir=results_dir){
  if (!is.null(v.effects)){
    # ***************
    # add gene id to data frame of variant effects and factorise the variant effect for each donor
    # ***************
    v.df = v.df[, as.character(unique(donor.list))]
    mutational.landscape = mutate(v.df, GeneID = rownames(v.df)) %>% select(GeneID, colnames(v.df))
    # ***************
    # order the Variant_Effects so that legend and colours displayed in required order 
    # ***************
    v.donor = v.effects[grep("D0", v.effects)]
    v.effects = c(v.effects[order(v.effects)])
    if (length(which(v.effects=="NA"))>0) {
      na.added=TRUE
      na.pos = which(v.effects=="NA")
      v.effects = v.effects[c(1:(na.pos-1),(na.pos+1):length(v.effects), na.pos)]
    } else {
      na.added = FALSE
    }
    for (i in 2:ncol(mutational.landscape)){
      NA.level = which(is.na(mutational.landscape[,i])==TRUE)
      if (length(NA.level)>0){
        mutational.landscape[c(NA.level),i] ="NA"
        if (!na.added){
          v.effects = c(v.effects, "NA")
          na.added = TRUE
        }
      }
      NA.char = which(mutational.landscape[,i]=="NA")
      if (length(NA.char)>0 && !na.added){
        v.effects = c(v.effects, "NA")
        na.added = TRUE
      }
      mutational.landscape[,i] = factor(mutational.landscape[,i], levels=v.effects)
    }
    # ***************
    # create blank and # Mutations in matrix column in format of melt function
    # ***************
    num.mutations = NULL
    for (i in 1:nrow(mutational.landscape)){
      if (rownames(mutational.landscape)[i]%in%c(top_annots)){ #"Status", 
        num.na = num.gain = num.loss = num.germ = 0
      } else {
        num.gain = length(which(mutational.landscape[i,]=="Gain"))
        num.loss = length(union(which(mutational.landscape[i,]== "cn2LOH"), which(mutational.landscape[i,]== "LOH")))                          
        num.na = length(which(mutational.landscape[i,]=="NA")) + length(which(is.na(mutational.landscape[i,])))
        num.germ = length(which(mutational.landscape[i,]=="germline"))
      }
      num.mutations = c(num.mutations, (ncol(mutational.landscape) -1 - num.na - num.gain - num.loss - num.germ))
    }
    names(num.mutations) = mutational.landscape$GeneID
    for (extra.gene in gene.list$Hugo_Symbol[grep("_", gene.list$Hugo_Symbol)]){
      num.mutations[extra.gene] = num.mutations[unlist(strsplit(extra.gene, "_"))[1]]
    }
    # ***************
    # write the data to be plotted to a txt file
    # ***************
    info.file = cbind(mutational.landscape, num.mutations)
    colnames(info.file)[ncol(info.file)] = "# Cases"
    if (!is.null(gene.order)){
      top.annots = rownames(mutational.landscape)[which(!rownames(mutational.landscape)%in%gene.list$Hugo_Symbol)]
      if (is.null(gene.order)){
        stop("You must supply a list of genes in the order you want plotted.")
      } else {
        info.file$GeneID = factor(info.file$GeneID, levels=c(top.annots, gene.order))
        info.file = info.file[order(info.file$GeneID),]
      }
    } else {
      info.file = info.file[order(info.file$`# Cases`, decreasing=TRUE),]
    }
    for (i in 1:ncol(info.file)){
      info.file[,i] = as.character(info.file[,i])
    }
    for (i in 1:nrow(info.file)){
      if (rownames(info.file)[i]%in%gene.list$Hugo_Symbol){
        for (k in 1:ncol(info.file)){
          if (colnames(info.file)[k]%in%donor.list){
            if (!is.na(info.file[i,k])){
              mutations = unlist(strsplit(as.character(info.file[i,k]), ';'))
              first = TRUE
              for (mut in mutations){
                mut = gsub("_2", "", gsub("_3", "", gsub("_4", "", mut)))
                if (length(grep(mut, SNV_mutations))>0){
                  aa = intersect(which(data$UniqSampleLabel==colnames(info.file)[k]),
                                 intersect(which(data$Hugo_Symbol==rownames(info.file)[i]),
                                           grep(mut, data$Plot_Class)))
                  if (length(aa)>1){
                    value = paste(data$Amino_Acid_Change[aa], collapse=";")
                  } else {
                    value = data$Amino_Acid_Change[aa]
                  }
                  if (first){
                    new_value = paste(mut, value, sep=";")
                  } else {
                    new_value = paste(new_value, paste(mut, value, sep=";"), sep=":")
                  }
                  first = FALSE
                } else {
                  if (first){
                    new_value = mut
                    first = FALSE
                  } else {
                    new_value = paste(new_value, mut, sep=":")
                  }
                }
              }
              info.file[i,k] = new_value
            }
          }
        }
      }
    }
    write_tsv(info.file, 
              file=paste(results.dir, v.filename, "_", ifelse(nchar(plot.type)>1, 
                                                              paste0(plot.type, "_cohort_", ncol(v.df), "_samples.tsv"), 
                                                              paste0(ncol(v.df), "_samples.tsv")), sep="")) 
    no.mutations = which(num.mutations==0)
    if (length(no.mutations)>0 && is.null(gene.order)){
      gene.list = gene.list[-which(gene.list$Hugo_Symbol%in%names(num.mutations[no.mutations])),]
      gene.list$Hugo_Symbol = droplevels(gene.list$Hugo_Symbol)
      num.mutations = num.mutations[-no.mutations]
      mutational.landscape = mutational.landscape[-no.mutations,]
    } else if (length(no.mutations)>0){
      num.mutations = num.mutations[-no.mutations]
      mutational.landscape = mutational.landscape[-no.mutations,]
    }
    mutational.landscape %<>% filter(GeneID%in%names(num.mutations))
    # ***************
    # use melt function to create a dataframe for the mutational landscape for plotting with ggplot
    # ***************
    mut.gene = reshape2::melt(mutational.landscape, id.vars = 'GeneID', variable.name = "DonorID")
    plot.data = mut.gene
    if (!is.null(gene.order)){
      if (is.null(gene.order)){
        stop("You must supply a list of genes in the order you want plotted.")
      } else {
        plot.data$GeneID = factor(plot.data$GeneID, 
                                  levels=c(top.annots, gene.order))
      }
    } else {
      plot.data$GeneID = factor(plot.data$GeneID, 
                                levels=c(as.character(names(num.mutations)[order(num.mutations, decreasing=TRUE)])))
    }
    v.cellularity = unique(plot.data[intersect(which(plot.data$GeneID=="Cellularity"), 
                                               which(plot.data$DonorID%in%colnames(v.df))),c("value")])
    if (length(v.cellularity)>0){
      v.effects = v.effects[grep(paste(v.cellularity, collapse="|"), v.effects, invert = TRUE)]
      v.cellularity = v.cellularity[order(v.cellularity)]
      v.cellularity = v.cellularity[union(grep("<50%", v.cellularity), union(grep("% to <", v.cellularity, fixed=TRUE), grep(">=90|Cell_NA",v.cellularity)))]
    }
    v.ploidy = unique(plot.data[intersect(which(plot.data$GeneID=="Ploidy"), 
                                          which(plot.data$DonorID%in%colnames(v.df))),c("value")])
    if (length(v.ploidy)>0){
      v.effects = v.effects[grep(paste(v.ploidy, collapse="|"), v.effects, invert = TRUE)]
      if (length(v.ploidy)>1){
        if (!"PLNA" %in% v.ploidy && !"PL>5"%in%v.ploidy){
          v.ploidy = v.ploidy[order(round(as.numeric(gsub("PL", "", v.ploidy))))]
        } else {
          v.ploidy.num = v.ploidy[-grep("PLNA|PL>5", v.ploidy)]
          v.ploidy.num = v.ploidy.num[order(round(as.numeric(gsub("PL", "", v.ploidy.num))))]
          if ("PL>5" %in% v.ploidy ){
            v.ploidy.num = c(v.ploidy.num, "PL>5")
          }
          if (("PLNA" %in% v.ploidy )){
            v.ploidy.num = c(v.ploidy.num, "PLNA")
          }
          v.ploidy = v.ploidy.num
        }
      }
    }
    if (status){
      if (length(which(is.na(unique(plot.data$value[which(plot.data$GeneID=="Status")]))))>0 | 
                 length(which(unique(plot.data$value[which(plot.data$GeneID=="Status")])=="Status_NA"))>0) {
        v.status = c(levels(study.data$Donor_class), "NA")
      } else {
        v.status = levels(study.data$Donor_class)
      }
      v.effects = v.effects[grep(paste(v.status, collapse="|"), v.effects, invert = TRUE)]
    } else {
      v.status = NULL
    }
    v.tmb = unique(plot.data[intersect(which(plot.data$GeneID=="TMB"), 
                                       which(plot.data$DonorID%in%colnames(v.df))),c("value")])
    if (length(v.tmb)>0){
      v.effects = v.effects[grep(paste(v.tmb, collapse="|"), v.effects, invert = TRUE)]
      v.tmb = v.tmb[order(v.tmb)]
      v.tmb = v.tmb[union(grep("<=5", v.tmb), union(grep(">5", v.tmb), union(grep(">10",v.tmb), union(grep(">15", v.tmb), grep("TMB_NA", v.tmb)))))]
    }
    v.msi = unique(plot.data[intersect(which(plot.data$GeneID=="MSI"), 
                                       which(plot.data$DonorID%in%colnames(v.df))),c("value")])
    if (length(v.msi)>0){
      v.effects = v.effects[grep(paste(v.msi, collapse="|"), v.effects, invert = TRUE)]
      v.msi = v.msi[order(v.msi)]
      v.msi = v.msi[union(grep("<3.5", v.msi), union(grep(">3.5", v.msi), grep("MSI_NA", v.msi)))]
    }
    v.meth.clust = unique(plot.data[intersect(which(plot.data$GeneID=="Meth_Cluster"), 
                                              which(plot.data$DonorID%in%colnames(v.df))),c("value")])
    if (length(v.meth.clust)>0){
      v.effects = v.effects[grep(paste(v.meth.clust, collapse="|"), v.effects, invert = TRUE)]
      if ("Meth_Cluster_NA" %in% v.meth.clust){
        v.meth.clust = v.meth.clust[-which(v.meth.clust=="Meth_Cluster_NA")]
        v.meth.clust = v.meth.clust[order(v.meth.clust)]
        v.meth.clust = c(v.meth.clust, "Meth_Cluster_NA")
      } else {
        v.meth.clust = v.meth.clust[order(v.meth.clust)]
      }
    }
    v.wgd = unique(plot.data[intersect(which(plot.data$GeneID=="WGD_Status"), 
                                       which(plot.data$DonorID%in%colnames(v.df))),c("value")])
    if (length(v.wgd)>0){
      v.effects = v.effects[grep(paste(v.wgd, collapse="|"), v.effects, invert = TRUE)]
      if ("WGD_NA" %in% v.wgd){
        v.wgd = v.wgd[-which(v.wgd=="WGD_NA")]
        v.wgd = v.wgd[order(v.wgd)]
        v.wgd = c(v.wgd, "WGD_NA")
      } else {
        v.wgd = v.wgd[order(v.wgd)]
      }
    }
    v.amp.arch = unique(plot.data[intersect(which(plot.data$GeneID=="Amp_Arch"), 
                                            which(plot.data$DonorID%in%colnames(v.df))),c("value")])
    if (length(v.amp.arch)>0){
      v.effects = v.effects[grep(paste(v.amp.arch, collapse="|"), v.effects, invert = TRUE)]
      if ("Amp_Arch_NA" %in% v.amp.arch){
        v.amp.arch = v.amp.arch[-which(v.amp.arch=="Amp_Arch_NA")]
        v.amp.arch = v.amp.arch[order(v.amp.arch)]
        v.amp.arch = c(v.amp.arch, "Amp_Arch_NA")
      } else {
        v.amp.arch = v.amp.arch[order(v.amp.arch)]
      }
    }
    v.treatment = unique(plot.data[intersect(which(plot.data$GeneID=="Treatment"), 
                                             which(plot.data$DonorID%in%colnames(v.df))),c("value")])
    if (length(v.treatment)>0){
      v.effects = v.effects[grep(paste(v.treatment, collapse="|"), v.effects, invert = TRUE)]
      if ("Treatment_NA" %in% v.treatment){
        v.treatment = v.treatment[-which(v.treatment=="Treatment_NA")]
        v.treatment = v.treatment[order(v.treatment)]
        v.treatment = c(v.treatment, "Treatment_NA")
      } else {
        v.treatment = v.treatment[order(v.treatment)]
      }
    }
    v.survival = unique(plot.data[intersect(which(plot.data$GeneID=="Survival"), 
                                            which(plot.data$DonorID%in%colnames(v.df))),c("value")])
    if (length(v.survival)>0){
      v.effects = v.effects[grep(paste(v.survival, collapse="|"), v.effects, invert = TRUE)]
      if ("Survival_NA" %in% v.survival){
        v.survival = v.survival[-which(v.survival=="Survival_NA")]
        v.survival = v.survival[order(v.survival)]
        v.survival = c(v.survival, "Survival_NA")
      } else {
        v.survival = v.survival[order(v.survival)]
      }
    }
    v.effects = c(v.effects, v.status, v.cellularity, v.ploidy, v.tmb, v.msi, v.meth.clust, v.wgd, v.amp.arch, v.treatment, v.survival)
    returned.data = list()
    returned.data[["Data"]]=plot.data
    returned.data[["Effects"]]=v.effects
    returned.data[["NumMutations"]] = num.mutations
    returned.data[["Cellularity"]] = v.cellularity
    returned.data[["Ploidy"]] = v.ploidy
    returned.data[["TMB"]] = v.tmb
    returned.data[["MSI"]] = v.msi
    returned.data[["Status"]] = v.status
    returned.data[["Meth_Cluster"]] = v.meth.clust
    returned.data[["WGD_Status"]] = v.wgd
    returned.data[["Amp_Arch"]] = v.amp.arch
    returned.data[["Treatment"]] = v.treatment
    returned.data[["Survival"]] = v.survival
    return(returned.data)
  } else {
    return(NULL)
  }
}

# ***************
# Function to set up colours to be used in the plots and to populate values
# ***************
get.plot.colours = function(plot.data=NULL, 
                            type, 
                            gene.list = Gene_list){
  if (is.null(plot.data)){
    stop("In order to define plot colours the plot data must not be NULL")
  }
  effect = plot.data$Effects
  num.mut = plot.data$NumMutations[which(names(plot.data$NumMutations)%in%gene.list$Hugo_Symbol)]
  num.annot = plot.data$NumMutations[which(!names(plot.data$NumMutations)%in%gene.list$Hugo_Symbol)]
  cellularity = plot.data$Cellularity
  ploidy = plot.data$Ploidy
  tmb = plot.data$TMB
  msi = plot.data$MSI
  status = plot.data$Status
  meth.clust = plot.data$Meth_Cluster
  wgd = plot.data$WGD_Status
  amp.arch = plot.data$Amp_Arch
  treat = plot.data$Treatment
  survival = plot.data$Survival
  # soc = plot.data$SOC
  # genetype = plot.data$GeneType
  blank.col = "#ffffff"
  na.col = "#f1f1f1"
  annot.na.col = "#f1f1f1"
  zero.col = "#dddddd"
  if (!is.null(effect)){
    # Set up colours for gains, losses, other variants and no value in the matrix, 
    blank.var = which(effect == "blank")
    just.na.var = which(effect == "NA")
    if (length(which(effect == "TotalLoss"))>0){ #&& type%in%c("CNV", "summary_all", "summary")
      just.totalloss.var = which(effect == "TotalLoss")
    } else {
      just.totalloss.var = NULL
    }
    # print(paste("just total loss var has ", length(just.totalloss.var), sep =""))
    if (length(which(effect == "HighGain"))>0){ #&& type%in%c("CNV", "summary_all", "summary"), union(which(effect == "HighGain(5-8)"), which(effect == "HighGain(>=9)"))
      just.highgain.var = which(effect == "HighGain")
    } else {
      just.highgain.var = NULL
    }
    # print(paste("just high gain var has ", length(just.highgain.var), sep=""))
    just.loss.var = switch(type,
                           summary_all = intersect(grep("Loss|LOH", effect), grep(";", effect, invert=TRUE)),
                           SNV = intersect(grep("_del", effect, ignore.case=TRUE), grep(";", effect, invert=TRUE)),
                           TOP = intersect(grep("Loss|LOH", effect), grep(";", effect, invert=TRUE)))
    just.gain.var = switch(type,
                           summary_all = union(which(effect == "FC<2"), which(effect == "FC>=2")),
                           SNV = intersect(grep("_ins", effect, ignore.case=TRUE), grep(";", effect, invert=TRUE)),
                           TOP = union(which(effect == "FC<2"), which(effect == "FC>=2")))
    num.var = which(effect %in% as.character(unique(num.mut)))
    num.annot.var = which(effect %in% as.character(unique(num.annot)))
    loss.var = switch(type,
                      summary_all = intersect(grep("Loss|LOH", effect), grep(";", effect)),
                      SNV = intersect(grep("_del", effect, ignore.case=TRUE), grep(";", effect)),
                      TOP = intersect(grep("Loss|LOH", effect), grep(";", effect)))
    gain.var = switch(type,
                      summary_all = union(intersect(grep("FC<", effect), grep(";", effect)), union(intersect(grep("FC>=2", effect), grep(";", effect)), intersect(grep("HighGain", effect), grep(";", effect)))),
                      SNV = intersect(grep("_ins", effect, ignore.case=TRUE), grep(";", effect)),
                      TOP = union(intersect(grep("HighGain", effect), grep(";", effect)), intersect(grep("FC", effect), grep(";", effect))))
    cell.var = which(effect %in% as.character(cellularity))
    ploid.var = which(effect %in% as.character(ploidy))
    tmb.var = which(effect %in% as.character(tmb))
    msi.var = which(effect %in% as.character(msi))
    meth.clust.var = which(effect %in% as.character(meth.clust))
    wgd.var = which(effect %in% as.character(wgd))
    amp.arch.var = which(effect %in% as.character(amp.arch))
    treat.var = which(effect %in% as.character(treat))
    survival.var = which(effect %in% as.character(survival))
    status.var = which(effect %in% as.character(status))
    var.set = unique(c(blank.var, just.na.var, just.totalloss.var, just.highgain.var, just.loss.var, just.gain.var, num.var,
                       num.annot.var, loss.var, gain.var, cell.var, ploid.var, tmb.var, status.var, 
                       msi.var, meth.clust.var, wgd.var, amp.arch.var, treat.var, survival.var))
    other.var = setdiff(seq(1:length(effect)), var.set)
    
    plot.colours = vector("character", length=length(effect))
    plot.colours[blank.var] = blank.col
    plot.colours[just.na.var] = na.col #plot.justs["NA"]
    # Assign Loss colours
    total.loss = as.numeric(sum(length(just.loss.var), length(just.totalloss.var)))
    # print(paste("total loss = ", total.loss, sep=""))
    loss.end = as.numeric(length(loss.var))
    plot.total.loss = colorRamp2(c(1, (total.loss+2)), c("#c4e8f8", "#0B5F83"), space="sRGB")
    plot.loss = colorRamp2(c(1, (loss.end+2)), c("#c3c9ee", "#0b1c89"), space="sRGB")
    # print(paste("plot loss has ", length(plot.loss(1:loss.end)), " colours", sep=""))
    # print(plot.loss(1:loss.end))
    just.totalloss.col = switch(type,
                                summary_all = plot.total.loss(total.loss+1),
                                SNV = NULL,
                                TOP = plot.total.loss(total.loss+1))
    # print(paste("just total loss colours has ", length(just.totalloss.col), " colours", sep=""))
    plot.colours[just.totalloss.var] = just.totalloss.col
    just.loss.col = switch(type,
                           summary_all = plot.total.loss(1:length(just.loss.var)),
                           SNV = plot.total.loss(1:length(just.loss.var)),
                           TOP = plot.total.loss(1:length(just.loss.var)))
    # print(paste("just loss colours has ", length(just.loss.col), " colours", sep=""))
    plot.colours[just.loss.var] = just.loss.col
    loss.col = switch(type,
                      summary_all = plot.loss(1:(length(loss.var))),
                      SNV = plot.loss(1:(length(loss.var))),
                      TOP = plot.loss(1:(length(loss.var))))
    # print(paste("loss colours has ", length(loss.col), " colours", sep=""))
    plot.colours[loss.var] = loss.col
    # Assign Gain colours
    total.gain = sum(length(just.gain.var), length(just.highgain.var))
    gain.end = ifelse(total.gain==1||total.gain==0, 3, total.gain+2)
    plot.total.gain = colorRamp2(c(1, (total.gain+2)), c("#FDC5C3", "#C80B04"), space="sRGB")
    plot.gain = colorRamp2(c(1, (length(gain.var)+2)), c("#FAD1DF", "#B31449"), space="sRGB")
    if (type %in% c("summary_all", "summary", "CNV")){
      if (length(just.highgain.var)>1){
        just.highgain.col = plot.total.gain((gain.end-1):gain.end)
      } else {
        just.highgain.col = plot.total.gain(gain.end)
      }
    } else {
      just.highgain.col = NULL
    }
    for (i in 1:length(just.highgain.var)){
      plot.colours[just.highgain.var[i]] = just.highgain.col[i]
    }
    just.gain.col = switch(type,
                           summary_all = plot.total.gain(2:(length(just.gain.var)+1)),
                           SNV = plot.total.gain(2:(length(just.gain.var)+1)),
                           TOP = plot.total.gain(2:(length(just.gain.var)+1)))
    plot.colours[just.gain.var] = just.gain.col
    gain.col = switch(type,
                      summary_all = plot.gain(1:(length(gain.var))),
                      SNV = plot.gain(1:(length(gain.var))),
                      TOP = plot.gain(1:(length(gain.var))))
    plot.colours[gain.var] = gain.col
    # Assign colours to other mutations
    other.names = effect[other.var]
    plot.other = c("darkseagreen2", "darkseagreen3", "gold1" , "gold2", "darkorchid3", "darkorchid4", "darkgoldenrod2",  "darkgoldenrod3", 
                   "magenta3", "magenta4", "navajowhite", "navajowhite3", "aquamarine3", "aquamarine4",  "indianred2", "indianred3", 
                   "bisque2", "bisque3", "tan2", "tan3", "#F781BF", "#C51B7D", "#7F00FF", "black")
    names(plot.other) = c("missense", "missense;+", "inframe", "inframe;+", "stop_gained", "stop_gained;+", "frameshift", "frameshift;+", 
                          "start_gained", "start_gained;+", "splice_site", "splice_site;+", "missense;SV", "SV", "deleterious", "deleterious;SV", 
                          "unknown", "unknown;SV", "Intragene_Fusion", "Intergene_Fusion", "other1", "other2", "other3", "multi-hit")
    num = 0
    for (o.var in other.var){
      if (effect[o.var]%in%names(plot.other)){
        plot.colours[o.var] = plot.other[effect[o.var]]
      } else {
        num = num + 1
        plot.colours[o.var] = plot.other[paste0("other", num)]
      }
    }
    # Assign number colours
    if (length(num.var)>1){
      plot.num = colorRamp2(c(0,max(num.mut)),c("#fee6e1","red3"), space="sRGB")
      for (i in 1:length(num.var)){
        plot.colours[num.var[i]] = plot.num(as.numeric(effect[num.var[i]]))
      }
    } else {
      plot.colours[num.var] = "#fee6e1"
    }
    
    plot.colours[num.annot.var] = "white"
    # Assign cellularity colours
    if (length(grep("TC", effect))>0){
      cell.colours = c()
      plot.cellularity = colorRamp2(c(0,100),c("#f7e4da","#793a1b"), space="sRGB") 
      for (i in 1:length(cellularity)){
        cell.colours = c(cell.colours, plot.cellularity(as.numeric(gsub("TC", "", cellularity[i]))))
      }
    } else {
      cell.colours =c("#dddddd", colorRamp2(seq(1,5, length=5),rev(brewer.pal(11, "BrBG")[1:5]))(2:5))
      names(cell.colours) = c("Cell_NA", "Cell_<50%", "Cell_50% to <75%", "Cell_75% to <90%", "Cell_>=90%")
      for (i in 1:length(cell.colours)){
        plot.colours[which(effect==names(cell.colours)[i])] = cell.colours[i]
      }
    }
    for (i in 1:length(cell.colours)){
      plot.colours[which(effect==names(cell.colours)[i])] = cell.colours[i]
    }
    # Assign ploidy colours
    if (length(grep("PL", effect))>0){
      if (length(grep("NA", ploidy))>0 && length(ploidy)==1){
        ploidy.colours = annot.na.col
      } else if (length(grep("NA", ploidy))>0){
        temp.ploidy = ploidy[-grep("NA", ploidy)]
        if (length(grep(">5", temp.ploidy))>0){
          plot.ploidy = colorRamp2(c(1,6), c("#e2d6f0", "#3a2159"), space="sRGB")
        } else {
          plot.ploidy = colorRamp2(c(1,max(as.numeric(gsub("PL", "", temp.ploidy)))), c("#e2d6f0", "#3a2159"), space="sRGB")
        }
        ploidy.colours = c()
        for (i in 1:length(temp.ploidy)){
          if (length(grep(">5", temp.ploidy[i]))>0){
            ploidy.colours = c(ploidy.colours, plot.ploidy(6))
          } else {
            ploidy.colours = c(ploidy.colours, plot.ploidy(as.numeric(gsub("PL", "", temp.ploidy[i]))))
          }
        }
        ploidy.colours = c(ploidy.colours, annot.na.col)
      } else {
        if (length(grep(">5", ploidy))>0){
          plot.ploidy = colorRamp2(c(1,6), c("#e2d6f0", "#3a2159"), space="sRGB")
        } else {
          plot.ploidy = colorRamp2(c(1,max(as.numeric(gsub("PL", "", ploidy)))), c("#e2d6f0", "#3a2159"), space="sRGB")
        }
        ploidy.colours = c()
        for (i in 1:length(ploidy)){
          if (length(grep(">5", ploidy[i]))>0){
            ploidy.colours = c(ploidy.colours, plot.ploidy(6))
          } else {
            ploidy.colours = c(ploidy.colours, plot.ploidy(as.numeric(gsub("PL", "", ploidy[i]))))
          }
        }
      }
      names(ploidy.colours) = ploidy
    } else {
      ploidy.colours = NULL
    }
    for (i in 1:length(ploidy.colours)){
      plot.colours[which(effect==names(ploidy.colours)[i])] = ploidy.colours[i]
    }
    # Assign tmb colours
    if ("TMB_NA" %in% tmb){
      if (length(tmb)>2){
        plot.tmb = colorRamp2(c(1,length(tmb)-1), c("#dff1fb", "#0c4564"), space="sRGB")
        tmb.colours = c(plot.tmb(1:(length(tmb)-1)), annot.na.col)
      } else if (length(tmb)==1) {
        tmb.colours = annot.na.col
      } else {
        tmb.colours = c("#0c4564", annot.na.col)
      }
    } else {
      plot.tmb = colorRamp2(c(1,length(tmb)), c("#dff1fb", "#0c4564"), space="sRGB")
      tmb.colours = c(plot.tmb(1:length(tmb)))
    }
    names(tmb.colours) = tmb
    for (i in 1:length(tmb.colours)){
      plot.colours[which(effect==names(tmb.colours)[i])] = tmb.colours[i]
    }
    # Assign msi colours
    plot.msi = colorRamp2(c(1,2), c("#dffbf5", "#0d6855"), space="sRGB")
    msi.colours = c(plot.msi(1:2), annot.na.col)
    names(msi.colours) = c("MSI_<3.5", "MSI_>3.5", "MSI_NA")
    for (i in 1:length(msi.colours)){
      plot.colours[which(effect==names(msi.colours)[i])] = msi.colours[i]
    }
    # Assign methylation cluster colours
    meth.clust.col = c("#49b05c", "#5108b4","#f58919", "#d9c80d", "#f188d5", "#610aa3", "#e62843", annot.na.col)
    names(meth.clust.col) = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "Meth_Cluster_NA")
    for (i in 1:length(meth.clust.col)){
      plot.colours[which(effect==names(meth.clust.col)[i])] = meth.clust.col[i]
    }
    # Assign wgd status colours
    wgd.col = c("lightgoldenrod1", "indianred1", annot.na.col)
    names(wgd.col) = c("HET", "DUP", "WGD_NA")
    for (i in 1:length(wgd.col)){
      plot.colours[which(effect==names(wgd.col)[i])] = wgd.col[i]
    }
    # Assign amplicon architecture status colours
    amp.arch.col = c("firebrick3", "cadetblue2", "khaki2", "coral1", "darkolivegreen2", "grey40", annot.na.col)
    names(amp.arch.col) = c("BFB", "Complex-non-cyclic", "ecDNA", "ecDNA;BFB", "ecDNA;Complex-non-cyclic", "None", "Amp_Arch_NA")
    for (i in 1:length(amp.arch.col)){
      plot.colours[which(effect==names(amp.arch.col)[i])] = amp.arch.col[i]
    }
    # Assign survival status colours
    survival.col = c("#8fd185", "#f45252", "khaki3", annot.na.col)
    names(survival.col) = c("Long", "Short", "Other", "Survival_NA")
    for (i in 1:length(survival.col)){
      plot.colours[which(effect==names(survival.col)[i])] = survival.col[i]
    }
    # Assign treatment colours
    treat.col = c("goldenrod1", "dodgerblue", annot.na.col)
    names(treat.col) = c("chemo/immuno", "other", "Treatment_NA")
    for (i in 1:length(treat.col)){
      plot.colours[which(effect==names(treat.col)[i])] = treat.col[i]
    }
    # Assign subtype colours
    status.levels = c("adeno", "NSCLC", "squamous", "sclc", "LC_neuro", "other_met", "other")
    level.cols = c("coral", "darkolivegreen1", "darkolivegreen4", "goldenrod1", 
                   "cadetblue3", "indianred1", "mediumorchid3", "darkseagreen1", "#555555", annot.na.col)[1:length(status.levels)] #, "darkseagreen3", "lightgoldenrod3" )
    level.cols = c(level.cols, annot.na.col)
    names(level.cols) = c(status.levels, "NA")
    status.colours = level.cols[status]
    for (i in 1:length(status.colours)){
      plot.colours[which(effect==names(status.colours)[i])] = status.colours[i]
    }
    names(plot.colours) = effect
    matrix.colours = list(plot.colours, just.gain.var, just.highgain.var, just.loss.var, 
                          just.totalloss.var, gain.var, loss.var, other.var, status.var, num.var, num.annot.var, 
                          cell.var, ploid.var, tmb.var, msi.var, meth.clust.var, wgd.var, amp.arch.var, treat.var, survival.var)
    names(matrix.colours) = c("plot.colours", "just.gain.var", "just.highgain.var", "just.loss.var", 
                              "just.totalloss.var", "gain.var", "loss.var", "other.var", "status.var", "num.var", "num.annot.var",
                              "cell.var", "ploid.var", "tmb.var", "msi.var", "meth.clust.var", "wgd.var", "amp.arch.var", "treat.var", "survival.var")
    return(matrix.colours)
  } else {
    return(NULL)
  }
}

# ***************
# Function to determine colours to use for the plot.
# ***************
get.value.colours = function(type,
                             plot.colours){
  switch(type,
         SNV = plot.colours[[type]]$plot.colours,
         TOP = plot.colours[[type]]$plot.colours,
         summary_all = plot.colours[[type]]$plot.colours)
}

# ***************
# function to create the plots
# ***************
create.matrix.plot = function(type,
                              plot.data,
                              plot.mut,
                              plot.colours,
                              plot.titles,
                              by.donor = FALSE,
                              study.data=NULL, 
                              top.annots=top_annots,
                              gene.list = NULL){
  # plot parameters
  bckgrnd.col="#ffffff"
  par(mai=c(0.1, 0.1, 0.1, 0.1))
  titlefontsize = 11
  lab.size.y = ifelse(length(levels(plot.data[[type]]$GeneID))<=40, 
                      9,
                      ifelse(length(levels(plot.data[[type]]$GeneID))<=60,
                             8, 5))
  
  lab.size.x = ifelse(length(levels(plot.data[[type]]$DonorID))<=40, 
                      8,
                      ifelse(length(levels(plot.data[[type]]$DonorID))<=80,
                             7, 8))
  if (report_gene=="MMR"){
    # overall TMB
    p.hist = filter(plot.data[[type]], GeneID=="Cellularity") %>% select(DonorID, value) %>% dplyr::rename(Cellularity=value)
    p.hist %<>% mutate(TMB = sapply(as.character(DonorID), function(x) study.data$Adj_TMB[grep(x, study.data$UniqueLabel)], USE.NAMES=FALSE),
                       Cell_col = sapply(Cellularity, function(x) plot.colours[which(names(plot.colours[[type]]$plot.colours)==x)],
                                         USE.NAMES=FALSE))
    p.hist$TMB[which(is.na(p.hist$TMB))] = 0
    p.hist$DonorID = factor(p.hist$DonorID, levels=unique(plot.data[[type]]$DonorID))
    p.annotations = filter(plot.data[[type]], GeneID%in%c(top.annots)) %>% filter(!GeneID=="TMB") #"Status", 
    p.annotations$GeneID = droplevels(p.annotations$GeneID)
    # add Indel TMB
    p.hist %<>% mutate(Indel_TMB = sapply(as.character(DonorID), function(x) study.data$Indel_TMB[grep(x, study.data$UniqueLabel)], USE.NAMES=FALSE))
    p.hist$Indel_TMB[which(is.na(p.hist$Indel_TMB))] = 0
  } else {
    p.annotations = filter(plot.data[[type]], GeneID%in%c(top.annots)) #"Status", 
    p.annotations$GeneID = droplevels(p.annotations$GeneID)
    p.hist = NULL
  }
  p.data = filter(plot.data[[type]], !GeneID%in%c(top.annots)) #"Status", 
  p.data$GeneID = droplevels(p.data$GeneID)
  if (type%in%c("summary_all", "TOP")){
    if ("Group"%in%colnames(gene.list)){
      p.data %<>% mutate(Group = sapply(as.character(p.data$GeneID), function(x) gene.list$Group[which(gene.list$Hugo_Symbol==x)]))
      p.data$Group = factor(p.data$Group, levels = unique(gene.list$Group))
      p.data$GeneID = as.character(p.data$GeneID)
      gene.levels = c()
      for (gp.level in levels(gene.list$Group)){
        gene.levels = c(gene.levels, gene.list$Hugo_Symbol[which(gene.list$Group==gp.level)])
      }
      gene.levels = gene.levels[which(gene.levels%in%unique(p.data$GeneID))]
      p.data$GeneID = factor(p.data$GeneID, levels = gene.levels)
    }
  }
  # create the different plots
  num.sample = length(unique(p.data$DonorID[grep(substr(p.data$DonorID[1], 1, 10), p.data$DonorID)]))
  vline.index = cumsum(table(study.data$Meth_Cluster))
  if (!is.null(gene.list) & "Group"%in%colnames(gene.list)){
    gene.list %<>% dplyr::arrange(desc(Group)) %>% dplyr::filter(Hugo_Symbol%in%unique(p.data$GeneID))
    gene.list$Group = factor(gene.list$Group, levels = unique(gene.list$Group))
    hline.index = cumsum(table(gene.list$Group))
  } else {
    hline.index = NULL
  }
  # gene plots
  if (!is.null(p.data) && nrow(p.data)>0){
    if ("Group" %in% colnames(p.data)){
      if (type=="summary_all") {
        p = list()
        p.gene.hist = list()
        for (gp in levels(p.data$Group)){
          pgrid.data = dplyr::filter(p.data, Group==gp)
          pgrid.data$GeneID = droplevels(pgrid.data$GeneID)
          plot.ghist.data = data.frame(GeneID=unique(pgrid.data$GeneID),
                                       NumMut = sapply(as.character(unique(pgrid.data$GeneID)), function(x) plot.mut[[type]][x]))
          plot.ghist.data$GeneID = factor(plot.ghist.data$GeneID, levels=rev(levels(pgrid.data$GeneID)))
          if (nrow(pgrid.data)>0){
            p[[gp]] = ggplot(pgrid.data, aes(DonorID, GeneID)) +
              geom_tile(aes(fill=value), colour="white") +
              scale_fill_manual(values=get.value.colours(type,
                                                         plot.colours), na.value=bckgrnd.col)
            y.label = ifelse(gp=="Cell_Cycle_Regulation", "Cell Cycle\nRegulation",
                             ifelse(gp=="Transcriptional_Regulation", "Transcriptional\nRegulation",
                                    ifelse(gp=="Receptor_Kinase_P13K_signalling", "Receptor Kinase\nP13K Signalling",
                                           ifelse(gp=="Notch_signalling_NE_differentiation", "NOTCH Signal\nNE Diff", 
                                                  ifelse(gp=="Other", "Other SCLC Genes", "Gene ID")))))
            p[[gp]] = p[[gp]] + xlab(ifelse(by.donor, "Donor Label", "Sample Label")) + ylab(y.label)
            p[[gp]] = p[[gp]] + theme(plot.title = element_blank())
            p[[gp]] = p[[gp]] + theme(axis.text.y=element_text(size=lab.size.y), 
                                      legend.position="none",
                                      panel.grid.major=element_blank(),
                                      panel.background=element_rect(fill=bckgrnd.col))
            p[[gp]] = p[[gp]] + ylim(rev(levels(pgrid.data$GeneID))) 
            if (y.label=="Gene ID"){
              if (!gp==levels(p.data$Group)[ceiling(length(levels(p.data$Group))/2)]){
                p[[gp]] = p[[gp]] + theme(axis.title.y=element_blank())
              }
            }
            for (ind in vline.index[1:(length(vline.index)-1)]){
              p[[gp]] = p[[gp]] + geom_segment(x=ind+0.5, xend = ind+0.5, y=(1-0.48),
                                               yend=length(levels(pgrid.data$GeneID))+0.48, color="black")
            }
            p.gene.hist[[gp]] = ggplot(plot.ghist.data, aes(x=GeneID, y=NumMut)) + 
              geom_bar(stat="identity", fill="lightblue") + ylim(0, max(plot.mut[[type]])) + coord_flip() +
              theme(axis.text.y=element_blank(), 
                    axis.title.y=element_blank(),
                    panel.grid.major.y=element_blank(),
                    panel.grid.minor.x=element_blank(),
                    panel.background = element_rect(fill="#f1f1f1"))
            
            if (!gp==levels(p.data$Group)[length(levels(p.data$Group))]){
              p[[gp]] = p[[gp]] + theme(axis.text.x=element_blank(),
                                        axis.title.x=element_blank(),
                                        axis.ticks=element_blank())
              p.gene.hist[[gp]] = p.gene.hist[[gp]] + theme(axis.text.x=element_blank(),
                                                            axis.title.x=element_blank(),
                                                            axis.ticks = element_blank())
            } else {
              p[[gp]] = p[[gp]] + theme(axis.text.x=element_text(angle=90, hjust=1, size=lab.size.x))
              p.gene.hist[[gp]] = p.gene.hist[[gp]] + theme(axis.ticks.y=element_blank())
            }
          } else {
            p[[gp]] = NULL
            p.gene.hist[[gp]] = NULL
          }
        }
        num.gene = length(unique(p.data$GeneID))
      } else if (type=="TOP") {
        p = list()
        p.gene.hist = list()
        for (gp in unique(p.data$Group)){
          pgrid.data = dplyr::filter(p.data, Group==gp)
          pgrid.data$GeneID = droplevels(pgrid.data$GeneID)
          plot.ghist.data = data.frame(GeneID=unique(pgrid.data$GeneID),
                                       NumMut = sapply(as.character(unique(pgrid.data$GeneID)), function(x) nrow(dplyr::filter(pgrid.data, GeneID==x))))
          plot.ghist.data$GeneID = factor(plot.ghist.data$GeneID, levels=rev(levels(pgrid.data$GeneID)))
          if (nrow(pgrid.data)>0){
            p[[gp]] = ggplot(pgrid.data, aes(DonorID, GeneID)) +
              geom_tile(aes(fill=value), colour="white") +
              scale_fill_manual(values=get.value.colours(type,
                                                         plot.colours), na.value=bckgrnd.col)
            y.label = gp
            p[[gp]] = p[[gp]] + xlab(ifelse(by.donor, "Donor Label", "Sample Label")) + ylab(y.label)
            p[[gp]] = p[[gp]] + theme(plot.title = element_blank())
            p[[gp]] = p[[gp]] + theme(axis.text.y=element_text(size=lab.size.y), 
                                      legend.position="none",
                                      panel.grid.major=element_blank(),
                                      panel.background=element_rect(fill=bckgrnd.col))
            p[[gp]] = p[[gp]] + ylim(rev(levels(pgrid.data$GeneID))) 
            if (y.label=="Gene ID"){
              if (!gp==levels(p.data$Group)[ceiling(length(levels(p.data$Group))/2)]){
                p[[gp]] = p[[gp]] + theme(axis.title.y=element_blank())
              }
            }
            for (ind in vline.index[1:(length(vline.index)-1)]){
              p[[gp]] = p[[gp]] + geom_segment(x=ind+0.5, xend = ind+0.5, y=(1-0.48),
                                               yend=length(levels(pgrid.data$GeneID))+0.48, color="black")
            }
            p.gene.hist[[gp]] = ggplot(plot.ghist.data, aes(x=GeneID, y=NumMut)) + 
              geom_bar(stat="identity", fill="lightblue") + ylim(0, max(sapply(as.character(unique(pgrid.data$GeneID)), function(x) nrow(dplyr::filter(pgrid.data, GeneID==x))))) + coord_flip() +
              theme(axis.text.y=element_blank(), 
                    axis.title.y=element_blank(),
                    panel.grid.major.y=element_blank(),
                    panel.grid.minor.x=element_blank(),
                    panel.background = element_rect(fill="#f1f1f1"))
            
            if (!gp==levels(p.data$Group)[length(levels(p.data$Group))]){
              p[[gp]] = p[[gp]] + theme(axis.text.x=element_blank(),
                                        axis.title.x=element_blank(),
                                        axis.ticks=element_blank())
              p.gene.hist[[gp]] = p.gene.hist[[gp]] + theme(axis.text.x=element_blank(),
                                                            axis.title.x=element_blank(),
                                                            axis.ticks = element_blank())
            } else {
              p[[gp]] = p[[gp]] + theme(axis.text.x=element_text(angle=90, hjust=1, size=lab.size.x))
              p.gene.hist[[gp]] = p.gene.hist[[gp]] + theme(axis.ticks.y=element_blank())
            }
          } else {
            p[[gp]] = NULL
            p.gene.hist[[gp]] = NULL
          }
        }
        num.gene = length(unique(p.data$GeneID))
        }
      } else {
        p = ggplot(p.data, aes(DonorID, GeneID)) +
          geom_tile(aes(fill=value), colour="white") +
          scale_fill_manual(values=get.value.colours(type,
                                                     plot.colours), na.value=bckgrnd.col)
        p = p + xlab(ifelse(by.donor, "Donor Label", "Sample Label")) + ylab("Gene ID")
        p = p + theme(plot.title = element_blank())
        p = p + theme(axis.text.x=element_text(angle=90, hjust=1, size=lab.size.x),
                      axis.text.y=element_text(size=lab.size.y), 
                      legend.position="none",
                      axis.ticks=element_blank(),
                      panel.grid.major=element_blank(),
                      panel.background=element_rect(fill=bckgrnd.col))
        p = p + ylim(rev(levels(p.data$GeneID)))
        for (ind in vline.index[1:(length(vline.index)-1)]){
          p = p + geom_segment(x=ind+0.5, xend = ind+0.5, y=(1-0.48),
                               yend=length(levels(p.data$GeneID))+0.48, color="black")
        }
        if (!is.null(hline.index)){
          for (ind in hline.index[1:(length(hline.index)-1)]){
            p = p + geom_segment(x=(1-0.48), xend = length(levels(p.data$DonorID))+0.48, 
                                 y=ind+0.5, yend=ind+0.5, color="black")
          }
        }
        num.gene = length(unique(p.data$GeneID))
        p.gene.hist = NULL
    }
  } else {
    p = NULL
    p.gene.hist = NULL
    num.gene = 0
  }
  # annotations plot
  if (!is.null(p.annotations) && nrow(p.annotations)>0){
    donor.col = plot.colours[[type]]$plot.colours[grep("D0", names(plot.colours[[type]]$plot.colours))]
    donor.pos = length(unique(p.annotations$GeneID))-1
    p.anno = ggplot(p.annotations, aes(DonorID, GeneID)) +
      geom_tile(aes(fill=value), colour="white") +
      scale_fill_manual(values=get.value.colours(type,
                                                 plot.colours), na.value=bckgrnd.col)
    p.anno = p.anno + ggtitle(paste(study_name, plot.titles[type], analysis_name, sep=" "))
    p.anno = p.anno + theme(plot.title = element_blank())
    p.anno = p.anno + theme(axis.text.x=element_blank(),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank(),
                            axis.text.y=element_text(size=lab.size.y), legend.position = "none",
                            axis.ticks = element_blank(),
                            panel.grid.major = element_blank(),
                            panel.background = element_rect(fill=bckgrnd.col))
    p.anno = p.anno + ylim(rev(levels(p.annotations$GeneID)))
    num = 1
    for (ind in vline.index[1:(length(vline.index)-1)]){
      p.anno = p.anno + geom_segment(x=ind+0.5, xend = ind+0.5, y=(1-0.48),
                           yend=length(levels(p.annotations$GeneID))+0.48, color="black")
    }
    num.anno = length(unique(p.annotations$GeneID))
    num.donor = length(unique(p.annotations$DonorID))
    for (i in 1:num.anno){
      p.anno = p.anno + geom_segment(x=0.5, xend = num.donor+0.5, y=i-0.5, 
                                     yend=i-0.5, color="white")
    }
  } else {
    p.anno = NULL
    num.anno = 0
  }
  # tmb hist
  if (!is.null(p.hist) && nrow(p.hist)>0){
    if ("TMB" %in% colnames(p.hist)){
      p.tmb = ggplot(p.hist, aes(x=DonorID, y=TMB)) + geom_bar(stat="identity", color=c("#527E96FF"), fill=c("#527E96FF")) 
      p.tmb = p.tmb + theme(axis.text.x=element_blank(),
                            axis.title.x=element_blank(),
                            axis.text.y=element_text(size=7), legend.position = "none",
                            axis.ticks = element_blank(),
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.y = element_line(size = 0.25),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_line(color="white", size=1),
                            panel.background = element_rect(fill = "#f1f1f1"))
    }
  } else {
    p.tmb = NULL
    p.tmb.ind = NULL
  }
  plot.space = ifelse(type=="summary_all", -0.17, -0.17)
  anno.height = ifelse(type=="summary_all", num.anno*1.2, ifelse(num.gene<60, 4, ifelse(num.gene<80, 3.5, ifelse(num.gene<120, 2.5, 1.25))))
  if (type%in%c("summary_all", "TOP") & !is.null(gene.list) & "Group"%in%colnames(gene.list)){
    gene.height = c()
    for (plot in names(p)) {
      plot.height = ifelse(!is.null(p[[plot]]), 
                           ifelse(!plot==names(p)[length(names(p))], 
                                  length(unique(p.data$GeneID[which(p.data$Group==plot)]))*1.2, 
                                  ifelse(type=="TOP", length(unique(p.data$GeneID[which(p.data$Group==plot)]))*1.2+2, length(unique(p.data$GeneID[which(p.data$Group==plot)]))*1.2+2.2)),
                           0)
      gene.height = c(gene.height, plot.height, plot.space)
    }
  } else {
    gene.height = ifelse(num.gene<60, 7, ifelse(num.gene<80, 7.5, ifelse(num.gene<120, 8, 9.5)))
  }
  if (!is.null(p.tmb) && !is.null(p.tmb.ind) && !is.null(p.anno) && !is.null(p)){
    plot.heights = c(4, plot.space, 4, plot.space, anno.height, plot.space, gene.height)
  } else if (!is.null(p) && !is.null(p.anno)){ 
    plot.heights = c(anno.height, plot.space, gene.height)
  } else if (!is.null(p.tmb) && !is.null(p.anno)){
    plot.heights = c(3, plot.space, anno.height)
  } else if (!is.null(p.tmb) && !is.null(p)){
    plot.heights = c(3, plot.space, gene.height)
  }
  p.list = list(tmb.plot = p.tmb, tmb.ind.plot = p.tmb.ind, anno.plot = p.anno, gene.plot = p, gene.hist = p.gene.hist, plot.heights = plot.heights)
  return(p.list)
}

# ***************
# function to create the legend
# ***************
create.legend = function(type,
                         plot.colours,
                         var.effects){
  plot.new()
  plot(NULL, xaxt='n', yaxt='n', bty='n', ylab='', xlab='', xlim=0:1, ylim=0:1)
  coord=par("usr")
  # set up the parameters for the legend
  legendtitlecex = 0.9 # 0.8
  legendtextcex = 0.7 # 0.6
  pointtype = 15
  pointcex = 1.5 # 1.2
  legend.ypos = c(0.97, 0.80, 0.56, 0.30)
  legend.xpos = c(-0.04, 0.25, 0.45, 0.59, 0.75, 0.88)
  legend.colours = get.value.colours(type,
                                     plot.colours)
  status = plot.colours[[type]]$status.var
  wgd = plot.colours[[type]]$wgd.var
  cellularity = plot.colours[[type]]$cell.var
  ploidy = plot.colours[[type]]$ploid.var
  tmb = plot.colours[[type]]$tmb.var
  msi = plot.colours[[type]]$msi.var
  meth.clust = plot.colours[[type]]$meth.clust.var
  other = plot.colours[[type]]$other.var
  cases = plot.colours[[type]]$num.var
  treatment = plot.colours[[type]]$treat.var
  survival = plot.colours[[type]]$survival.var
  amp.arch = plot.colours[[type]]$amp.arch.var
  losses = switch(type,
                  summary_all = union(plot.colours[[type]]$just.loss.var, 
                                      union(plot.colours[[type]]$just.totalloss.var, plot.colours[[type]]$loss.var)),
                  SNV = union(plot.colours[[type]]$just.loss.var, plot.colours[[type]]$loss.var),
                  TOP = union(plot.colours[[type]]$just.loss.var, 
                              union(plot.colours[[type]]$just.totalloss.var, plot.colours[[type]]$loss.var)))
  if (type=="TOP"){
    if (length(losses)>0){
      losses = c(losses, which(var.effects[[type]]=="NA"))
    }
  }
  
  gains = switch(type,
                 summary_all = union(plot.colours[[type]]$just.highgain.var, union(plot.colours[[type]]$just.gain.var, plot.colours[[type]]$gain.var)),
                 SNV = union(plot.colours[[type]]$just.gain.var, plot.colours[[type]]$gain.var),
                 TOP = union(plot.colours[[type]]$just.highgain.var, union(plot.colours[[type]]$just.gain.var, plot.colours[[type]]$gain.var)))
  
  # start legend
  num.pos = 1
  if (type=="TOP" & length(c(losses, gains))>0){
    legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[1], legend=c(""),
           title =ifelse(length(gains)>0, "LOH & Gains", "LOH"),
           col=c(legend.colours[gains], legend.colours[losses]),
           cex=legendtitlecex, bty="n")
  } else if (!type=="TOP" & length(c(losses, gains))>0) {
    legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[1], legend=c(""),
           title ="Copy Number Variants",
           col=c(legend.colours[gains], legend.colours[losses]),
           cex=legendtitlecex, bty="n")
  }
  if (type=="TOP") {
    var.effects[[type]][losses] = gsub("LOH", "YES", gsub("NA", "NO", var.effects[[type]][losses]))
  }
  if (length(gains)>0 || length(losses)>0){
    legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[1],
           legend=c(gsub("_Mutation", "", var.effects[[type]][gains]), 
                   gsub("_Mutation", "", var.effects[[type]][losses])),
           pch=pointtype, pt.cex=pointcex, title ="",
           col=c(legend.colours[gains], legend.colours[losses]),
           cex=legendtextcex, bty="n")
  }
  num.pos = num.pos + 1
  legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[1], legend=c(""),
           title ="Variant Classification",
           col=c(legend.colours[gains], legend.colours[losses]),
           cex=legendtitlecex, bty="n")
  if (length(other)>0){
    if (length(gains)>0 || length(losses)>0){
      legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[1],
             legend=c(gsub("_Mutation", "", var.effects[[type]][other])),
             pch=pointtype, pt.cex=pointcex, title ="",
             col=legend.colours[other],
             cex=legendtextcex, bty="n")
    } else {
      legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[1],
             legend=c(gsub("_Mutation", "", var.effects[[type]][other])),
             pch=pointtype, pt.cex=pointcex, title ="",
             col=legend.colours[other],
             cex=legendtextcex, bty="n")
    }
  }
  # Status, Methylation cluster and Survival 
  leg.lengths = unlist(lapply(list(status, meth.clust, survival), function(x) length(x)))
  names(leg.lengths) = c("Status", "Meth Cluster", "Survival")
  if (length(which(leg.lengths>0))>0){
    num.pos = num.pos + 1
    num.y.pos = 1
    for (i in 1:length(leg.lengths)){
      if (leg.lengths[i]>0){
        if (names(leg.lengths[i])=="Status"){
          leg.var = status
        } else if (names(leg.lengths[i])=="Meth Cluster"){
          leg.var = meth.clust
        } else if (names(leg.lengths[i])=="Survival"){
          leg.var = survival
        } 
        legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[num.y.pos], legend=c(""),
               title =names(leg.lengths[i]),
               col=c(legend.colours[leg.var]),
               cex=legendtitlecex, bty="n")
        legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[num.y.pos],
               legend=c(gsub("Cell_", "", gsub("T", "Group", var.effects[[type]][leg.var], fixed = TRUE))),
               pch=pointtype, pt.cex=pointcex, title="",
               col=c(legend.colours[leg.var]),
               cex=legendtextcex, bty="n")
        num.y.pos = num.y.pos + 1
      }
    }
  }
  # Treatment, Survival, Ploidy & MSI 
  leg.lengths = unlist(lapply(list(treatment, cellularity, ploidy, msi), function(x) length(x)))
  names(leg.lengths) = c("Treatment", "Cellularity", "Ploidy", "MSI")
  if (length(which(leg.lengths>0))>0){
    num.pos = num.pos + 1
    num.y.pos = 1
    for (i in 1:length(leg.lengths)){
      if (leg.lengths[i]>0){
        if (names(leg.lengths[i])=="Treatment"){
          leg.var = treatment
        } else if (names(leg.lengths[i])=="Cellularity"){
          leg.var = cellularity
        } else if (names(leg.lengths[i])=="Ploidy"){
          leg.var = ploidy
        } else if (names(leg.lengths[i])=="MSI"){
          leg.var = msi
        } 
        legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[num.y.pos], legend=c(""),
               title = ifelse(names(leg.lengths[i])=="Cellularity", "Tumour Content", names(leg.lengths[i])),
               col=c(legend.colours[leg.var]),
               cex=legendtitlecex, bty="n")
        legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[num.y.pos],
               legend=c(gsub("MSI_", "", gsub("PL", "", gsub("Cell_", "", var.effects[[type]][leg.var])))),
               pch=pointtype, pt.cex=pointcex, title="",
               col=c(legend.colours[leg.var]),
               cex=legendtextcex, bty="n")
        num.y.pos = num.y.pos + 1
      }
    }
  }
  # Whole genome duplication status, TMB and Amp Architecture
  leg.lengths = unlist(lapply(list(wgd, tmb, amp.arch), function(x) length(x)))
  names(leg.lengths) = c("WGD Status", "TMB", "Amp Architecture")
  if (length(which(leg.lengths>0))>0){
    num.pos = num.pos + 1
    num.y.pos = 1
    for (i in 1:length(leg.lengths)){
      if (leg.lengths[i]>0){
        if (names(leg.lengths[i])=="TMB"){
          leg.var = tmb
        } else if (names(leg.lengths[i])=="WGD Status"){
          leg.var = wgd
        } else if (names(leg.lengths[i])=="Amp Architecture"){
          leg.var = amp.arch
        }
        legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[num.y.pos], legend=c(""),
               title =names(leg.lengths[i]),
               col=c(legend.colours[leg.var]),
               cex=legendtitlecex, bty="n")
        legend(x=coord[2]*legend.xpos[num.pos], y=coord[4]*legend.ypos[num.y.pos],
               legend=c(gsub("TMB_", "", gsub("WGD_", "", gsub("Amp_Arch_", "", var.effects[[type]][leg.var])))),
               pch=pointtype, pt.cex=pointcex, title="",
               col=c(legend.colours[leg.var]),
               cex=legendtextcex, bty="n")
      }
      num.y.pos = num.y.pos + 1
    }
  }
}
