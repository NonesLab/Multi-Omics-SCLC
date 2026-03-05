# ***************
# This script is a repository of all the functions
# required when running the Mutations signature contribution script
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
read.in.data = function(pair.read.file,
                        full.study=TRUE, 
                        WGS.study=FALSE,
                        WES.study=FALSE){
  study.data = as.data.frame(fread(pair.read.file))
  
  # for cromwell limit entries in files that have an analysis
  if (full.study){
    study.data %<>% dplyr::filter(!is.na(analysis)) %>% dplyr::filter(!analysis=="")
  } else {
    study.data %<>% dplyr::filter(!is.na(analysis)) %>% dplyr::filter(!analysis=="")
    if (WGS_only){
      study.data %<>% dplyr::filter((is.na(testCaptureKit)|testCaptureKit=="")) %>% dplyr::filter((is.na(controlCaptureKit)|controlCaptureKit==""))
    } else if (WES_only) {
      study.data %<>% dplyr::filter(!is.na(testCaptureKit)) %>% dplyr::filter(!testCaptureKit=="") %>% dplyr::filter(!is.na(controlCaptureKit)) %>% dplyr::filter(!controlCaptureKit=="")
    }
  }
  # replace any spaces with underscores 
  # (this is a fix for the one donor record from one study that doesn't comply)
  # Donor Label
  study.data$donorLabel = gsub(" ", "_", study.data$donorLabel)
  # Sample Label
  study.data$testSampleLabel = gsub(" ", "_", study.data$testSampleLabel)
  study.data$controlSampleLabel = gsub(" ", "_", study.data$controlSampleLabel)
  data.sets = list(study=study.data)
  return(data.sets)
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
  get.ending = function(sample.vector,
                        ending.type){
    pos = ifelse(ending.type == "seq",
                 which(colnames(study.data)=="testCaptureKit"),
                 ifelse(ending.type == "platform",
                        which(colnames(study.data)=="testSequencingPlatform"),
                        which(colnames(study.data)=="donorLabel")))
    # message(paste("The position of comparison is:", pos))
    ending = ifelse(ending.type == "seq", 
                    ifelse(length(grep("custom", sample.vector[pos], ignore.case=TRUE))>0,
                           "Custom",
                           ifelse(length(grep("methyl", sample.vector[pos], ignore.case=TRUE))>0,
                                  "Methyl",
                                  ifelse(length(grep("panel|amplicon", sample.vector[pos], ignore.case=TRUE))>0,
                                         "Panel",
                                         ifelse(length(grep("exo", sample.vector[pos], ignore.case=TRUE))>0,
                                                "WES", "WGS")))),
                    ifelse(ending.type=="platform",
                           ifelse(as.character(sample.vector[pos])%in%c("NovaSeq", "HiSeq2000", "HiSeq2500", "HiSeqXTen", "MiSeq", 
                                                                        "NextSeq", "HiSeqXFive", "HiSeq4000", "HiSeq1500", "HiSeq3000"), 
                                  "Illumina",
                                  ifelse(as.character(sample.vector[pos])%in%c("BGIseq500", "MGIseq2000"), "BGI", "wrong")), 
                           ""))
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
                                                                                  "_", get.ending(x, ending.type), sep="")))
  } else {
    study.data %<>% mutate(UniqueLabel = apply(study.data, 1, function(x) paste(x[which(colnames(study.data)=="testSampleLabel")], 
                                                                                  "_", get.ending(x, ending.type), sep="")))
  }
  study.data$UniqueLabel[grep("dilute", study.data$testSampleLabel)] = paste(study.data$donorLabel[grep("dilute", study.data$testSampleLabel)], "WES_D", sep="_")
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
  # Sample Label
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
    sample.data = as.data.frame(fread(paste(samp.dir.name, samp.filename, sep=""), stringsAsFactors=FALSE))
  } else if (report && is.null(samp.dir.name) || is.null(samp.filename)) {
    stop("If you wish to report subgroups you must supply a directory and filename where the information can be found.")
  }
  if (!is.null(sample.data)){
    samp.col = names(which(apply(sample.data, 2, function(x) length(which(x %in% study.data[,samp.id]))>0)==TRUE))
    subgroups = data.frame(subgroup = sample.data[, subgp.var],
                           PID = sample.data[,samp.col])
    subgroup.df = data.frame(PID = vector(mode="character"),
                             subgroup = vector(mode="character"))
    for (i in 1:nrow(subgroups)){
      sample.labels = study.data$UniqueLabel[grep(subgroups$PID[i], study.data$UniqueLabel)]
      if (length(sample.labels)==1){
        subgroup.df = rbind(subgroup.df, data.frame(PID=sample.labels, subgroup=subgroups$subgroup[i]))
      } else {
        for (label in sample.labels){
          subgroup.df = rbind(subgroup.df, data.frame(PID=label, subgroup=subgroups$subgroup[i]))
        }
      }
    }
    if (length(which(is.na(subgroup.df$subgroup))>0)){
      message("All NAs found in the input subgroups were replaced with 'Other' ")
      subgroup.df$subgroup[which(is.na(subgroup.df$subgroup))] = "Other"
    }
    if (subgp.var=="Cat"){
      subgroup.df$subgroup = factor(subgroup.df$subgroup, levels=c("adeno", "non_scc", "squamous", "sclc", "LC_neuro", "other_met"))
      level_cols = c("coral", "darkolivegreen1", "darkolivegreen4", "lightgoldenrod1", 
                     "cadetblue", "indianred1", "mediumorchid3", "darkseagreen1", "#555555") 
    } else if (subgp.var=="Meth_Cluster"){
      subgroup.df$subgroup = factor(subgroup.df$subgroup, levels=c("T1", "T2", "T3", "T4"))
      level_cols = c("#49b05c", "#5108b4","#f58919", "#d9c80d", "#f188d5", "#610aa3", "#e62843")
    }
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
                              samp.id){
  samp.study.data = create.subgroup.data(report = report,
                                         samp.dir.name = samp.dir.name,
                                         samp.filename = samp.filename,
                                         subgp.var = subgp.var,
                                         study.data = study.data,
                                         samp.id = samp.id)
  sample.data = as.data.frame(fread(paste(samp.dir.name, samp.filename, sep=""), stringsAsFactors=FALSE))
  colnames(sample.data)[which(colnames(sample.data)==samp.id)] = "PID"
  sample.data$PID = sapply(sample.data$PID, function(x) paste(x, "WGS", sep="_"))
  sample.data = merge(samp.study.data, sample.data, by="PID")
  return(sample.data)
}

# *************
# function to read in the nominated analysis file and create sample data
# *************
read.in.study.data = function(pair.read=analysis_file,
                              full.study=full_study,
                              WGS.study=WGS_only,
                              WES.study=WES_only,
                              report = report_subgroups,
                              samp.dir.name=sample_data_dir,
                              samp.filename=sample_data_filename,
                              subgp.var=subgp_colname,
                              samp.id = samp_id){
  data.sets = read.in.data(pair.read.file=pair.read,
                           full.study=full.study, 
                           WGS.study=WGS.study,
                           WES.study=WES.study)
  study.data = data.sets$study
  study.data$donorLabel = gsub("-", "_", study.data$donorLabel)
  # *************
  # Filter summary files to samples of interest
  # *************
  samp.of.interest = as.data.frame(fread(paste(samp.dir.name, samp.filename, sep="")))
  study.data %<>% filter(donorLabel%in%samp.of.interest[,samp.id], testSampleLabel%in%samp.of.interest$testSampleLabel)
  
  # *************
  # creating an unique sample identifier for reporting
  # *************
  study.data = create.unique.label(study.data=study.data)
  # remove diluted samples for 100 Illumina exome study
  pattern = paste(paste(paste("WES", seq(2,10,1), sep="_"), collapse="|"), paste(paste("WGS", seq(2,10,1), sep="_"), collapse="|"), sep="|")
  study.data %<>% filter(!grepl(pattern, UniqueLabel))
  study.data$UniqueLabel = gsub("WES_1", "WES", gsub("WGS_1", "WGS", study.data$UniqueLabel))
  study.data %<>% filter(!grepl("D01_19_004_WG", UniqueLabel))
  
  # *************
  # create the subgroups for the analyses
  # *************
  sample.data = create.sample.data(report = report,
                                   samp.dir.name=samp.dir.name,
                                   samp.filename=samp.filename,
                                   subgp.var=subgp.var,
                                   study.data=study.data,
                                   samp.id = samp.id)
  study.data = list(study=study.data, sample=sample.data)
}


# *************
# Function to assign signature matrix and colours
# *************
get.sig.matrix.details = function(sbs.matrix.name,
                                  indel.matrix.name,
                                  sbs.cutoff,
                                  indel.cutoff,
                                  work.dir = paste0(study_dir, "Signatures/COSMIC_Info/"),
                                  package = NULL){
  if (package == "YAPSA"){
    # *************
    # load YAPSA data sets
    # *************
    data(sigs) 
    data(cutoffs)
    data(sigs_pcawg)
    data(cutoffs_pcawg)
  }
  
  # *************
  # Declare customised colours for SBS and Indel signatures
  # *************
  UV = c(brewer.pal(n=9, name="Pastel1")[6],
         brewer.pal(n=8, name="Pastel2")[6],
         brewer.pal(n=8, name="Accent")[4],
         brewer.pal(n=12, name="Set3")[12],
         brewer.pal(n=8, name="Set2")[6])
  # lime green - 2
  clockwise = c(brewer.pal(n=8, name="Set2")[5],
                brewer.pal(n=12, name="Set3")[7])
  # reds - 2
  apobec = c(brewer.pal(n=9, name="Reds")[6],
             brewer.pal(n=9, name="YlOrRd")[6])
  # blues - 5
  mmr = c(brewer.pal(n=12, name="Set3")[5],
          brewer.pal(n=9, name="PuBuGn")[6],
          brewer.pal(n=12, name="Paired")[1],
          brewer.pal(n=11, name="RdYlBu")[9],
          brewer.pal(n=11, name="RdBu")[8])
  # hot pinks - 5
  pols = c(brewer.pal(n=9, name="Set1")[8],
           brewer.pal(n=9, name="RdPu")[c(4,5)],
           brewer.pal(n=8, name="Set2")[4],
           brewer.pal(n=8, name="Dark2")[4])
  # purples - 2
  mmr_pols = c(brewer.pal(n=12, name="Set3")[10],
               brewer.pal(n=9, name="Set1")[4])
  # Sig.3 orange, Sig.18 navy, Sig.30 & Sig.36 dark blue - 4
  other = c(brewer.pal(n=9, name="YlOrRd")[5],
            brewer.pal(n=9, name="YlGnBu")[7],
            brewer.pal(n=11, name="PiYG")[1],
            brewer.pal(n=9, name="GnBu")[9])
  # orange/pinks - 3
  tobacco = c(brewer.pal(n=12, name="Paired")[5],
              brewer.pal(n=12, name="Set3")[4],
              brewer.pal(n=9, name="Pastel1")[1])
  # greens - 7
  treatment = c(brewer.pal(n=11, name="RdYlGn")[c(9,10,11)],
                brewer.pal(n=12, name="Paired")[4],
                brewer.pal(n=11, name="PRGn")[c(9,10,11)]
                )
  # blue greens - 7
  exposures = c(brewer.pal(n=11, name="BrBG")[c(7,8,9,10)],
                brewer.pal(n=11, name="Spectral")[9],
                brewer.pal(n=8, name="Pastel2")[1],
                brewer.pal(n=8, name="Dark2")[1])
  # - 18
  unknown = c(brewer.pal(n=8, name="Set2")[7],
              brewer.pal(n=8, name="Pastel2")[c(2,7)],
              brewer.pal(n=9, name="Pastel1")[c(5,7)],
              brewer.pal(n=8, name="Accent")[c(3,7)],
              brewer.pal(n=12, name="Paired")[c(7,12)],
              brewer.pal(n=9, name="Set1")[7],
              brewer.pal(n=8, name="Dark2")[c(6,7)],
              brewer.pal(n=11, name="BrBG")[c(1,2,3,4,5)],
              "#ad6a1a")
  unvalidated = c("gray20", "gray23", "gray26", "gray29", "gray32", "gray35", "gray38", "gray41", "gray44", 
                  "gray47", "gray50", "gray53", "gray56", "gray59", "gray62", "gray65", "gray68", "gray71", "gray91")
  # SBS colours
  all.sig.colours = c(clockwise[1], apobec[1], other[1], tobacco[1], clockwise[2], mmr[1],
                      UV[1], UV[1], UV[2], UV[3], UV[4], unknown[1], pols[1], 
                      pols[2], pols[2], pols[3], pols[4], pols[5], treatment[1], unknown[2],
                      apobec[2], mmr_pols[1], mmr[2], unknown[3], unknown[4], unknown[4], unknown[5],
                      other[2], unknown[6], mmr_pols[2], mmr[3], exposures[1], unknown[7],
                      exposures[2], treatment[2], mmr[4], unvalidated[1], unknown[8], tobacco[2],
                      other[3], treatment[3], treatment[4], unknown[9], unknown[10], treatment[5],
                      other[4], unknown[11], UV[5], unknown[12], unknown[13], unknown[14],
                      exposures[3], unvalidated[2], mmr[5], unvalidated[3], unvalidated[4], unvalidated[5],
                      unvalidated[6], unvalidated[7], unvalidated[8], unvalidated[9], unvalidated[10], unvalidated[11],
                      unvalidated[12], unvalidated[13], unvalidated[14], unvalidated[15], unvalidated[16], unvalidated[17],
                      unvalidated[18], exposures[4], exposures[5], treatment[6], treatment[7], exposures[6],
                      unknown[15], exposures[7], unknown[16], tobacco[3], unknown[17], unknown[18], unvalidated[19], "#d9f2ed")
  names(all.sig.colours) = c("Sig.1", "Sig.2", "Sig.3", "Sig.4", "Sig.5", "Sig.6", 
                             "Sig.7", "Sig.7a", "Sig.7b", "Sig.7c", "Sig.7d", "Sig.8", 
                             "Sig.9", "Sig.10", "Sig.10a", "Sig.10b", "Sig.10c", "Sig.10d", 
                             "Sig.11", "Sig.12", "Sig.13", "Sig.14", "Sig.15", "Sig.16", 
                             "Sig.17", "Sig.17a", "Sig.17b", "Sig.18", "Sig.19", "Sig.20", 
                             "Sig.21", "Sig.22", "Sig.23", "Sig.24", "Sig.25", "Sig.26", 
                             "Sig.27", "Sig.28", "Sig.29", "Sig.30", "Sig.31", "Sig.32", 
                             "Sig.33", "Sig.34", "Sig.35", "Sig.36", "Sig.37", "Sig.38", 
                             "Sig.39", "Sig.40", "Sig.41", "Sig.42", "Sig.43", "Sig.44", 
                             "Sig.45", "Sig.46", "Sig.47", "Sig.48", "Sig.49", "Sig.50", 
                             "Sig.51", "Sig.52", "Sig.53", "Sig.54", "Sig.55", "Sig.56", 
                             "Sig.57", "Sig.58", "Sig.59", "Sig.60", "Sig.84", "Sig.85", 
                             "Sig.86", "Sig.87", "Sig.88", "Sig.89", "Sig.90", "Sig.91", 
                             "Sig.92", "Sig.93", "Sig.94", "unassigned", "LT30")
  # Indel colours
  indel.colours = c(brewer.pal(n=12, name="Paired")[c(3,4)], tobacco[2], unknown[1], unknown[2], other[1], 
                    mmr[1], mmr_pols[1], unknown[6], unknown[7], unknown[8], unknown[9],
                    UV[5], unknown[10], unknown[15], unknown[16], mmr_pols[2], exposures[3], unvalidated[19], "#d9f2ed")
  indel.names = c("ID1", "ID2", "ID3", "ID4", "ID5", "ID6", "ID7",  
                  "ID8", "ID9", "ID10", "ID11", "ID12", "ID13", "ID14",
                  "ID15", "ID16", "ID17", "ID18", "unassigned", "LT30")
  names(indel.colours) = indel.names
  
  comb.colours = c(clockwise[1], apobec[1], other[1], tobacco[1], mmr[1], UV[5], pols[1], pols[5], mmr_pols[1], mmr_pols[2],
                   unknown[1], brewer.pal(n=9, name="Blues")[9], brewer.pal(n=9, name="Blues")[9], unknown[7], exposures[4], treatment[1], treatment[2], treatment[3], treatment[4], treatment[5],
                   treatment[6], exposures[1], treatment[7], exposures[2], exposures[7], unvalidated[1], unvalidated[18], "#d9f2ed")
  names(comb.colours) = c("Age", "APOBEC", "HR_Def", "Tobacco", "MMR", "UV", "POLE", "POLD1", "POLE_MMR", "POLD1_MMR",
                          "ROS_unclear", "ROS_MUTYH", "MUTYH", "NTHL1", "AID", "Temozolomide", "Aristolochic", "Aflatoxin", "Chemo", "Platinum", 
                          "Azathioprine", "Haloalkane", "Thiopurine", "Colibactin", "Duocarmycin", "Artifact", "Unknown", "LT30")
  # Assign matrix of SBS signatures to parameter 'sig.matrix'
  if (package=="YAPSA"){
    if (sbs.matrix.name=="COSMIC v2 (30 Sig)"){
      # import SBS signatures provided by YAPSA
      sig.matrix = AlexCosmicValid_sig_df
      colnames(sig.matrix) = gsub("AC", "Sig.", colnames(sig.matrix))
      # read in the cosmic signature optimised cutoffs or assign cutoff to 'sig_cutoff'
      if (!is.numeric(sbs.cutoff)){
        sig.cutoff = cutoffCosmicValid_abs_df[6,]
        names(sig.cutoff) = gsub("AC", "Sig.", names(sig.cutoff))
      } else {
        sig.cutoff = sbs.cutoff
      }
      # Assign df of signature characteristics to parameter 'sig_ind_df'
      sig.ind.df = AlexCosmicValid_sigInd_df
      sig.ind.df$sig = gsub("AC", "Sig.", sig.ind.df$sig)
      } else if (sbs.matrix.name=="COSMIC v3 Real"){
        # choose matrix
        sig.matrix = PCAWG_SP_SBS_sigs_Real_df
        colnames(sig.matrix) = gsub("SBS", "Sig.", colnames(sig.matrix)) 
        # read in the cosmic signature optimised cutoffs or assign cutoff to 'sig_cutoff'
        if (!is.numeric(sbs.cutoff)){
          sig.cutoff = cutoffPCAWG_SBS_WGSWES_realPid_df[10,]
          names(sig.cutoff) = gsub("SBS", "Sig.", names(sig.cutoff))
        } else {
          sig.cutoff = sbs.cutoff
        }
        # Assign df of signature characteristics to parameter 'sig_ind_df'
        sig.ind.df = PCAWG_SP_SBS_sigInd_Real_df
        sig.ind.df$sig = gsub("SBS", "Sig.", sig.ind.df$sig)
      } else if (sbs.matrix.name=="COSMIC v3 Artif"){
        # choose matrix
        sig.matrix = PCAWG_SP_SBS_sigs_Artif_df
        sig.matrix = sig.matrix[,grep("SBS84|SBS85", colnames(sig.matrix), invert = TRUE)]
        colnames(sig.matrix) = gsub("SBS", "Sig.", colnames(sig.matrix))
        # read in the cosmic signature optimised cutoffs or assign cutoff to 'sig_cutoff'
       if (!is.numeric(sbs.cutoff)){
          sig.cutoff = cutoffPCAWG_SBS_WGSWES_artifPid_df[16,]
          names(sig.cutoff) = gsub("SBS", "Sig.", names(sig.cutoff))
        } else {
          sig.cutoff = sbs.cutoff
        }
        # Assign df of signature characteristics to parameter 'sig_ind_df'
        sig.ind.df = PCAWG_SP_SBS_sigInd_Artif_df
        sig.ind.df = sig.ind.df[grep("SBS84|SBS85", sig.ind.df$sig, invert = TRUE),]
        sig.ind.df$sig = gsub("SBS", "Sig.", sig.ind.df$sig)
      } else {
        # Add custom matrix with cutoffs and colours here
        new.matrix = read.delim(paste(work.dir, "COSMIC_v3.2_SBS_GRCh37.txt", sep=""),
                                sep="\t", header=TRUE)
        sig.matrix = as.matrix(new.matrix[,2:ncol(new.matrix)], nrow=nrow(new.matrix))
        rownames(sig.matrix) = new.matrix$Type
        rownames(sig.matrix) = paste(substr(rownames(sig.matrix),3,5), " ", 
                                     substr(rownames(sig.matrix),1,1), 
                                     substr(rownames(sig.matrix),3,3), 
                                     substr(rownames(sig.matrix),7,7),sep="")
        sig.matrix = sig.matrix[order(rownames(sig.matrix)),]
        colnames(sig.matrix) = gsub("SBS", "Sig.", colnames(sig.matrix)) 
        sig.matrix = as.data.frame(sig.matrix)
        if (!is.numeric(sbs.cutoff)){
          stop("There are no optimised cutoffs for this signature matrix.")
        } else {
          sig.cutoff = sbs.cutoff
        }
        # Assign df of signature characteristics to parameter 'sig_ind_df'
        sig.ind.df = read.delim(paste(work.dir, "COSMIC_v3.2_SBS_ind_df.tsv", sep=""),
                                sep="\t", header=TRUE)
        sig.ind.df$sig = gsub("SBS", "Sig.", sig.ind.df$sig)
      }
  } else if (package=="MutationalPatterns"){
    sig.cutoff = sbs.cutoff
    new.matrix = read.delim(paste(work.dir, "COSMIC_v3.2_SBS_GRCh37.txt", sep=""),
                            sep="\t", header=TRUE)
    sig.matrix = as.matrix(new.matrix[,2:ncol(new.matrix)], nrow=nrow(new.matrix))
    rownames(sig.matrix) = new.matrix$Type
    rownames(sig.matrix) = paste(substr(rownames(sig.matrix),3,5), " ",
                                 substr(rownames(sig.matrix),1,1),
                                 substr(rownames(sig.matrix),3,3),
                                 substr(rownames(sig.matrix),7,7),sep="")
    sig.matrix = sig.matrix[order(rownames(sig.matrix)),]
    colnames(sig.matrix) = gsub("SBS", "Sig.", colnames(sig.matrix))
  }
  # *************
  # Assign customised colours and signature order for reporting 
  # *************
  sig.order = colnames(sig.matrix)
  sig.colours = all.sig.colours[c(colnames(sig.matrix), "unassigned", "LT30")]
  names(sig.colours) = c(colnames(sig.matrix), "unassigned", "LT30")
  if (package=="YAPSA"){
    sig.ind.df$colour = sig.colours[match(sig.ind.df$sig, names(sig.colours))]
  } else {
    sig.ind.df = NULL
  }
  
  # Assign matrix and cutoffs for indel signature analysis
  if (package=="YAPSA") {
    if (indel.matrix.name=="PCAWG Indel v3 (16 Sig)"){
      indel.sig.matrix = PCAWG_SP_ID_sigs_df
      if (!is.numeric(indel.cutoff)){
        indel.sig.cutoff = cutoffPCAWG_ID_WGS_Pid_df[3,]
      } else {
        indel.sig.cutoff = indel.cutoff
      }
      indel.sig.ind.df = PCAWG_SP_ID_sigInd_df
    } else {
      # Add custom matrix with cutoffs and colours here
      new.matrix = read.delim(paste(work.dir, "COSMIC_v3.2_ID_GRCh37.txt", sep=""),
                              sep="\t", header=TRUE)
      indel.sig.matrix = as.matrix(new.matrix[,2:ncol(new.matrix)], nrow=nrow(new.matrix))
      rownames(indel.sig.matrix) = new.matrix$Type
      indel.sig.matrix = as.data.frame(indel.sig.matrix)
      if (!is.numeric(indel.cutoff)){
        stop("Have only implemented optimum cutoffs for PCAWG Indel v3 (16 Sig) matrix supplied with YAPSA.")
      } else {
        indel.sig.cutoff = indel.cutoff
      }
      indel.sig.ind.df = read.delim(paste(work.dir, "COSMIC_v3.2_ID_ind_df.tsv", sep=""),
                                    sep="\t", header=TRUE)
    }
  } else if (package=="MutationalPatterns"){
    indel.sig.matrix = get_known_signatures(muttype="indel")
    indel.sig.cutoff = indel.cutoff
  }
  # *************
  # Assign customised colours for reporting 
  # *************
  if (package == "YAPSA"){
    indel.sig.ind.df$colour = indel.colours[match(indel.sig.ind.df$sig, names(indel.colours))]
  } else {
    indel.sig.ind.df = NULL
  }
  return(list(sbs_matrix = sig.matrix,
              cutoff = sig.cutoff,
              sbs_order = sig.order,
              sbs_colour = sig.colours,
              sbs_sig_df = sig.ind.df,
              indel_matrix = indel.sig.matrix,
              indel_cutoff = indel.sig.cutoff,
              indel_colour = indel.colours,
              indel_sig_df = indel.sig.ind.df,
              combined_cols = comb.colours))

}


# ******************
# Function to copy over vcfs into a working folder
# ******************
copy.vcf.data = function(study.data=study_data,
                         file.end=vcf_file_ending,
                         work.dir=working_dir, 
                         temp.dir=temp_directory){
  # populate the directory
  for (i in 1:nrow(study.data)){
    # read in the data
    if (!is.na(study.data$analysisPath[i])){
      vcf.dir = paste(work.dir, "vcfs/", sep="")
      if(!dir.exists(vcf.dir)){
        dir.create(vcf.dir, recursive=TRUE)
      }
      message(paste("Copying across vcf for:", study.data$donorLabel[i]))
      message(paste("File located in:", study.data$analysisPath[i]))
      copy.file = list.files(study.data$analysisPath[i], paste(file.end, "$",sep=""), full.names=TRUE)
      if (length(copy.file)==1){
        system(paste("cp ", copy.file, " ", vcf.dir, study_data$UniqueLabel[i], ".vcf.gz", sep=""))
      }
    }
  }
}

# ******************
# Function to build SNP type and transition for study
# ******************
build.vcf.data = function(study.data=study_data,
                          SNV.file.end=SNV_HC_file,
                          temp.dir=temp_directory){
  vcf.data = data.frame(Sample=vector(mode="character"),
                        chr=vector(mode="character"),
                        pos=vector(mode="character"),
                        ref=vector(mode="integer"),
                        alt=vector(mode="numeric"),
                        stringsAsFactors = FALSE
  )
  
  # populate the data
  for (i in 1:nrow(study.data)){
    # read in the data
    if (!is.na(study.data$analysisPath[i])){
      if (length(dir(study.data$analysisPath[i], pattern=paste(SNV.file.end, "$", sep="")))==1){
        SNV.sample.file = paste(study.data$analysisPath[i], "/", dir(study.data$analysisPath[i], pattern=paste(SNV.file.end, "$", sep="")), sep="")
        message(paste("Reading in vcf for:", study.data$donorLabel[i]))
        message(paste("File located in:", study.data$analysisPath[i]))
        system(paste("zcat ", SNV.sample.file, " |sed '/^##/d' > ", 
                     temp.dir, "snv.tsv", sep=""), intern=TRUE)
        type.data = read.delim(paste(temp.dir, "snv.tsv", sep=""), check.names=FALSE,
                               header=TRUE, sep="\t", na.strings = "", stringsAsFactors=FALSE)
        system(paste("rm ", temp.dir, "snv.tsv", sep=""))
      }
    }
    colnames(type.data)[which(colnames(type.data)=="#CHROM")] = "Chromosome"
    type.data = filter(type.data, grepl("chr", Chromosome))
    type.data %<>% mutate(Sample=study.data$UniqueLabel[i],
                          chr=type.data$Chromosome,
                          pos=type.data$POS,
                          ref=type.data$REF,
                          alt=type.data$ALT) %>% select(Sample,
                                                        chr,
                                                        pos,
                                                        ref,
                                                        alt,
                                                        FILTER,
                                                        INFO,
                                                        FORMAT)
    vcf.data = rbind(vcf.data, type.data)
  }
  return(list(vcf=vcf.data))
}

# ****************
# Function to find create the signature data
# ****************
get.signatures = function(sample.vcf.input,
                          count.mat=NULL,
                          WES.flag = WES_only,
                          WGS.flag = WGS_only,
                          ref.genome=ref_genome,
                          sig.matrix,
                          packages=c("YAPSA"),
                          analysis.type = c("SBS", "Indel"),
                          sig.cutoff=sig_cutoff,
                          sig.ind.df = NULL){
  if (!("YAPSA"%in%packages)){
    stop("function built specifically for package parameter 'YAPSA'")
  } 
  if (missing(sig.matrix)){
    stop("You must supply a matrix of known signatures for the supervised analysis.")
  }
  if (!(analysis.type=="SBS" | analysis.type=="Indel")){
    stop("analysis.type parameter should be one of 'SBS', or 'Indel'")
  } 
  if (WES.flag && WGS.flag){
    stop("This analysis can only run on all WGS data or all WES data.  Please set the flag in the declarations.")
  }
  if (!WES.flag && !WGS.flag){
    stop("This function needs to know whether your data is all WGS samples or all WES samples  Please set the flag in the declarations.")
  }
  flag = ifelse(WGS.flag, "WGS", "WES")
  sig.list = list()
  for (package in packages){
    if (package=="YAPSA"){
      sample.input = sample.vcf.input
      if (!is.data.frame(sample.input)) {
        stop(paste("for", package, "package the sample.input parameter should be a data.frame of vcf values."))
      }
      if (package=="YAPSA"){
        if (is.null(ref.genome)){
          stop(paste("For", package, "package you must supply a reference genome."))
        }
        if (is.null(sig.ind.df)){
          message("No signature information 'sig.ind.df' has been supplied so you will not be able to run confidence levels on the results.")
        }
        if (analysis.type=="SBS"){
          colnames(sample.input) = c("PID", "CHROM", "POS", "REF", "ALT")
          sample.input$change = attribute_nucleotide_exchanges(sample.input)
          sample.input %<>% arrange(PID, CHROM, POS)
          sample.input = annotate_intermut_dist_cohort(sample.input, in_PID.field="PID")
        } 
        message("Have annotated YAPSA data.")
      } else {
        stop(paste("could not initialise data for", package))
      }
    }
    message("Starting analysis and creating Signatures")
    if (package=="YAPSA"){
      # create data and extract signatures
      if (analysis.type=="SBS"){
        if (is.null(count.mat)){
          sig.list[[package]] = create_mutation_catalogue_from_df(this_df=sample.input,
                                                                  this_seqnames.field="CHROM", 
                                                                  this_start.field="POS", 
                                                                  this_end.field="POS", 
                                                                  this_PID.field="PID", 
                                                                  this_refGenome=ref.genome, 
                                                                  this_wordLength=3,
                                                                  this_verbose=FALSE)
          if (flag=="WES"){
            data("targetCapture_cor_factors")
            target.capture = "Agilent6withoutUTRs"
            cor.list = targetCapture_cor_factors[[target.capture]]
            sig.list[[package]]$matrix_old = sig.list[[package]]$matrix
            sig.list[[package]]$matrix = normalizeMotifs_otherRownames(sig.list[[package]]$matrix_old,
                                                                       cor.list$rel_cor)
          }
        } else {
          sig.list[[package]]$matrix = count.mat
        }
      } else if (analysis.type=="Indel"){
        if (is.null(count.mat)){
          sig.list[[package]]$matrix = create_indel_mutation_catalogue_from_df(in_dat = sample.input,
                                                                               in_signature_df = sig.matrix,
                                                                               in_verbose = FALSE)
        } else {
          sig.list[[package]]$matrix = count.mat
        }
      } else {
        stop("analysis.type should be either 'SBS', or 'Indel'")
      }
      if (analysis.type=="SBS"){
        yapsa.cutoff = sig.cutoff[[package]]
      } else if (analysis.type=="Indel"){
        yapsa.cutoff = sig.cutoff
      }
      sig.list[[package]]$LCD_cutoff = LCD_complex_cutoff_combined(in_mutation_catalogue_df=as.data.frame(sig.list[[package]]$matrix),
                                                                   in_signatures_df=sig.matrix,
                                                                   in_cutoff_vector=as.vector(yapsa.cutoff),
                                                                   in_sig_ind_df = sig.ind.df)
    }
  }
  return(sig.list)
} 

# ****************
# Function to find remove samples with no extracted signatures
# ****************
filter.yapsa.result = function(sig.list){
  lose.samples = names(which(colSums(sig.list$exposures)==0))
  if (length(lose.samples)>0){
    keep.samples = names(which(colSums(sig.list$exposures)>0))
    sig.list$exposures = sig.list$exposures[,keep.samples]
    sig.list$norm_exposures = sig.list$norm_exposures[,keep.samples]
    sig.list$residual_catalogue = sig.list$residual_catalogue[,keep.samples]
    sig.list$confidence %<>% filter(sample%in%keep.samples)
  }
  return(sig.list)
}
  
  
# ****************
# Function to find create the confidence levels of the extracted signatures
# ****************
get.confidence.levels = function(sig.list,
                                 package=c("YAPSA")){
  if (!(package=="YAPSA")){
    stop("This function is for the output from the 'YAPSA' R package only.")
  }
  sig.list[[package]]$LCD_cutoff$perPID$confidence = variateExp(in_catalogue_df = as.data.frame(sig.list[[package]]$matrix),
                                                                in_sig_df = sig.list[[package]]$LCD_cutoff$perPID$signatures,
                                                                in_exposures_df = sig.list[[package]]$LCD_cutoff$perPID$exposures
                                                                )
  sig.list[[package]]$LCD_cutoff$cohort$confidence = variateExp(in_catalogue_df = as.data.frame(sig.list[[package]]$matrix),
                                                                in_sig_df = sig.list[[package]]$LCD_cutoff$cohort$signatures,
                                                                in_exposures_df = sig.list[[package]]$LCD_cutoff$cohort$exposures
                                                                )
  sig.list[[package]]$LCD_cutoff$consensus$confidence = variateExp(in_catalogue_df = as.data.frame(sig.list[[package]]$matrix),
                                                                   in_sig_df = sig.list[[package]]$LCD_cutoff$consensus$signatures,
                                                                   in_exposures_df = sig.list[[package]]$LCD_cutoff$consensus$exposures
                                                                   )
  return(sig.list)
}

# ******************
# Function to plot a summary of the contribution of signatures
# ******************
plot.contrib.summary = function (contribution, 
                                 signatures, 
                                 index=c(), 
                                 subtypes.df=NULL,
                                 coord.flip=FALSE, 
                                 plot.title="",
                                 mode="relative", 
                                 palette=c(),
                                 summary.plot = FALSE){
  if (!(mode == "relative" | mode == "absolute")) 
    stop("mode parameter should be either 'relative' or 'absolute'")
  if (length(index > 0)) {
    contribution = contribution[, index]
  }
  Sample = NULL
  Contribution = NULL
  Signature = NULL
  if (mode == "relative") {
    m_contribution = reshape2::melt(contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
  } else {
    if (missing(signatures)) 
      stop(paste("For contribution plotting in mode 'absolute':", 
                 "also provide signatures matrix"))
    contribution = contribution[,which(colSums(contribution)>0)]
    sig.order = colnames(signatures)[which(colnames(signatures)%in%rownames(contribution))]#sig_order[which(sig_order%in%rownames(contribution))]
    sig.order = sig.order[grep("unassigned", sig.order, invert = TRUE)]
    signatures = matrix(as.matrix(signatures[,c(sig.order)]), ncol=length(sig.order), dimnames=list(rownames(signatures), sig.order))
    total_signatures = colSums(signatures)
    if ("unassigned" %in% rownames(contribution)){
      total_signatures = c(total_signatures, unassigned=1)
    }
    abs_contribution = contribution * total_signatures
    m_contribution = reshape2::melt(abs_contribution)
    colnames(m_contribution) = c("Signature", "Sample", "Contribution")
  }
  if (!is.null(subtypes.df)) {
    if (!is.data.frame(subtypes.df))
        stop("Subtypes should be a data frame of sample id 'PID' and sample class 'subgroup' and colour 'col'.")
    m_contribution$Subgroup = NA
    m_contribution$colours = NA
    for (sample in unique(m_contribution$Sample)){
      subgroup = as.character(subtypes.df$subgroup[which(subtypes.df$PID==sample)])
      subgroup_col = subtypes.df$col[which(subtypes.df$PID==sample)]
      m_contribution$Subgroup[which(m_contribution$Sample==sample)] = subgroup
      m_contribution$colours[which(m_contribution$Sample==sample)] = subgroup_col
    }
    if (is.null(levels(subtypes.df$subgroup))){
      m_contribution$Subgroup = factor(m_contribution$Subgroup, levels=unique(m_contribution$Subgroup))
    } else {
      m_contribution$Subgroup = factor(m_contribution$Subgroup, levels=levels(subtypes.df$subgroup))
    }
  }
  if (mode == "relative") {
    plot = ggplot(m_contribution, 
                  aes(x=factor(Sample), 
                      y=Contribution, 
                      fill=factor(Signature), 
                      color=factor(Signature),
                      order=Sample)) 
    if (summary.plot){
      plot = plot + geom_bar(stat="identity", width=0.8) #+ ylim(0,101)
    } else {
      plot = plot + geom_bar(position="fill", stat="identity",colour="white", width=0.98) 
    }
    plot = plot + labs(x="", y="Relative contribution") 
    if (!summary.plot){
      plot = plot + theme_bw()
    }
  } else {
    plot = ggplot(m_contribution, 
                  aes(x=factor(Sample), 
                      y=Contribution, 
                      fill=factor(Signature), 
                      color=factor(Signature),
                      order=Sample)) + 
      geom_bar(stat="identity", colour="white", width=0.98) + 
      labs(x="", y="Absolute contribution") +  theme_bw() 
  }
  
  plot = plot + ggtitle(plot.title)
  # determine number of column in legend.
  num.sig = length(unique(m_contribution$Signature))
  num.col = ifelse(num.sig<=25, 1, ifelse(num.sig<=50, 2, ifelse(num.sig<=75, 3, ifelse(num.sig<=100, 4 ,ifelse(num.sig<=75, 5, 6)))))
  # plot = plot + guides(fill=guide_legend(ncol=num.col), col=guide_legend(ncol=num.col))
  if (length(palette) > 0){
    plot = plot + scale_fill_manual(name="Signature", values=palette[as.character(unique(m_contribution$Signature))])
    plot = plot + scale_color_manual(name="Signature", values=palette[as.character(unique(m_contribution$Signature))])
  } else {
    plot = plot + scale_fill_discrete(name="Signature")
    plot = plot + scale_color_discrete(name="Signature")
  }
  plot = plot + theme(panel.grid.minor.x=element_blank(), 
                      panel.grid.major.x=element_blank(),
                      panel.background = element_blank(),
                      plot.title=element_text(size=title_size2),
                      axis.title.y=element_text(size=ylab_size),
                      axis.text.x=element_text(angle=60, hjust=1, size=xsize3))
  if (summary.plot){
    plot = plot + theme(axis.text.x=element_blank(),
                        axis.title.x=element_blank(),
                        axis.ticks=element_blank(),
                        legend.title=element_text(size=10),
                        legend.key.size=unit(legend_key_size2,"mm"),
                        legend.text=element_text(size=8))
  } else {
    plot = plot + theme(legend.title=element_text(size=legend_title_size2),
                        legend.key.size=unit(legend_key_size3,"mm"),
                        legend.text=element_text(size=legend_txt_size3))
  }
  if (!is.null(subtypes.df)){
    y_pos= ifelse(mode=="relative", -0.025, -1*(min(colSums(contribution))*0.075))
    y_pos_lab= ifelse(mode=="relative", 1.025, (max(colSums(contribution))*1.025))
    x_start = 0.5
    line_size=2
    if (!summary.plot){
      subtypes.df$subgroup = droplevels(subtypes.df$subgroup)
      for (level in levels(subtypes.df$subgroup)){
        x_end = x_start + length(unique(subtypes.df$PID[which(subtypes.df$subgroup==level)]))
        plot = plot + geom_segment(x=x_start, xend = x_end, y=y_pos, 
                                   yend=y_pos, color=unique(subtypes.df$col[which(subtypes.df$subgroup==level)]), linewidth=line_size)
        plot = plot + geom_label(label=level, x=(x_start + x_end)/2, y=y_pos_lab, size=3,
                                 label.padding = unit(0.1, "lines"), label.size=0.1, colour="black", 
                                 fill=unique(subtypes.df$col[which(subtypes.df$subgroup==level)]))
        x_start = x_end
      }
    } else {
      for (test.sample in colnames(contribution)){
        # level = as.character(unique(subtypes.df$PID[which(subtypes.df$PID==test.sample)]))
        x_end = x_start + 1
        plot = plot + geom_segment(x=x_start, xend=x_end, y=y_pos, yend=y_pos, color=subtypes.df$col[which(subtypes.df$PID==test.sample)], linewidth=line_size)
        x_start = x_end
      }
    }
  }
  if (coord.flip){
    plot = plot + coord_flip() + xlim(rev(levels(factor(m_contribution$Sample))))
  } 
  return(plot)
}

# ******************
# Function to plot a summary of the confidence levels of signatures
# ******************
myplotExposuresConfidence = function (in_complete_df, 
                                      in_subgroups_df, 
                                      in_sigInd_df,
                                      plot.title="") {
  sig_colour_vector <- c("slategrey", in_sigInd_df$colour)
  names(sig_colour_vector) <- c("total", as.character(in_sigInd_df$sig))
  in_complete_df$sample <- factor(in_complete_df$sample, levels = in_subgroups_df$PID[order(in_subgroups_df$index)])
  in_complete_df$subgroup <- in_subgroups_df$subgroup[match(in_complete_df$sample, 
                                                            in_subgroups_df$PID)]
  exposure_plot <- in_complete_df[which(in_complete_df$sig != "total"), ] %>% 
    ggplot(aes(x = sample, 
               y = exposure, 
               fill = sig)) + 
    geom_bar(stat = "identity") + 
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) + 
    facet_grid(sig ~ subgroup, scales = "free_x", space = "free_x", switch = "x") + 
    theme_grey() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=xsize2),
                         axis.text.y = element_text(size=xsize2),
                         panel.border = element_rect(fill = NA, colour = "black"), 
                         strip.background = element_rect(colour = "black"), 
                         legend.position = "none") + scale_fill_manual(values = sig_colour_vector)
  subgroup_aggregate_df <- aggregate(col ~ subgroup, data = in_subgroups_df, 
                                     FUN = head, 1)
  exposure_g <- ggplot_gtable(ggplot_build(exposure_plot))
  stripr <- which(grepl("strip-b", exposure_g$layout$name))
  for (i in stripr) {
    fill_index <- which(grepl("rect", exposure_g$grobs[[i]]$grobs[[1]]$childrenOrder))
    name_index <- which(grepl("text", exposure_g$grobs[[i]]$grobs[[1]]$childrenOrder))
    text_index <- which(grepl("text", exposure_g$grobs[[i]]$grobs[[1]]$children[[name_index]]))
    subgroup_name = exposure_g$grobs[[i]]$grobs[[1]]$children[[name_index]]$children[[text_index]]$label
    exposure_g$grobs[[i]]$grobs[[1]]$children[[fill_index]]$gp$fill <- 
      subgroup_aggregate_df$col[which(subgroup_aggregate_df$subgroup==subgroup_name)]
  }
  samp.order = c()
  textb = which(grepl("axis-b", exposure_g$layout$name))
  for (text_index in textb){
    title_index = which(grepl("title", exposure_g$grobs[[text_index]]$children$axis$grobs))
    label_index = which(grepl("text", exposure_g$grobs[[text_index]]$children$axis$grobs[[title_index]]$children))
    labels = exposure_g$grobs[[text_index]]$children$axis$grobs[[title_index]]$children[[label_index]]$label
    samp.order = c(samp.order, labels)
  }
  total_p <- in_complete_df[which(in_complete_df$sig == "total"),] %>% 
    ggplot(aes(x = sample, y = exposure, fill = sig)) + 
    geom_bar(stat = "identity", fill = "slategrey") + 
    ggtitle(plot.title) + 
    facet_grid(sig ~ subgroup, scales = "free_x", space = "free_x", switch = "x") + 
    theme_grey() + theme(axis.text.x = element_blank(), axis.ticks = element_blank(), 
                         axis.text.y = element_text(size=xsize2),
                         axis.title.x = element_blank(), axis.title.y = element_blank(), 
                         plot.title = element_text(size=title_size2),
                         panel.border = element_rect(fill = NA, colour = "black"), 
                         strip.background = element_rect(colour = "black"), strip.text.x = element_blank(), 
                         legend.position = "none")
  total_g <- ggplotGrob(total_p)
  all_g <- rbind(total_g, exposure_g)
  return(list(grob_all = all_g, order=samp.order))
}

