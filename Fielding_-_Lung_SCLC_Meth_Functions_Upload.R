# ***************
# This script is a repository of all the functions
# required when running the DNA Methylation Probe Analysis
# ***************

# ***************
# function to  read in the study summary spreadsheets for the mutations analysis summary
# ***************
read.in.data = function(analysis.file){
  study.data = as.data.frame(fread(analysis.file))
  # limit entries in files that have an analysis
  study.data %<>% dplyr::filter(!is.na(analysis)) %>% dplyr::filter(!analysis=="")
  # replace any spaces with underscores 
  study.data$donorLabel = gsub(" ", "_", study.data$donorLabel)
  # Sample Label
  study.data$testSampleLabel = gsub(" ", "_", study.data$testSampleLabel)
  study.data$controlSampleLabel = gsub(" ", "_", study.data$controlSampleLabel)
  data.set = list(study=study.data)
  return(data.set)
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
  BGI.seq = c("BGIseq500", "MGIseq2000")
  Illumina.seq = c("NovaSeq", "HiSeq2000", "HiSeq2500", "HiSeqXTen", "MiSeq", 
                   "NextSeq", "HiSeqXFive", "HiSeq4000", "HiSeq1500", "HiSeq3000")
  get.ending = function(sample.vector,
                        ending.type){
    pos = ifelse(ending.type == "seq" || ending.type == "all",
                 which(colnames(study.data)=="testCaptureKit"),
                 ifelse(ending.type == "platform",
                        which(colnames(study.data)=="testSequencingPlatform"),
                        which(colnames(study.data)=="donorLabel")))
    # message(paste("The position of comparison is:", pos))
    ending = ifelse(ending.type == "seq" || ending.type == "all", 
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
                                  ifelse(as.character(sample.vector[pos])%in%BGI.seq, "BGI", "wrong")), 
                           ""))
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
    if (length(unique(study.data$donorLabel))==nrow(study.data)){
      study.data %<>% mutate(UniqueLabel = donorLabel)
    } else {
      study.data %<>% mutate(UniqueLabel = apply(study.data, 1, function(x) paste0(x[which(colnames(study.data)=="donorLabel")], 
                                                                                   ifelse(nchar(get.ending(x, ending.type))>0, "_", ""), 
                                                                                   get.ending(x, ending.type))))
    }
  } else {
    if (length(unique(study.data$testSampleLabel))==nrow(study.data)){
      study.data %<>% mutate(UniqueLabel = testSampleLabel)
    } else {
      study.data %<>% mutate(UniqueLabel = apply(study.data, 1, function(x) paste0(x[which(colnames(study.data)=="testSampleLabel")], 
                                                                                   ifelse(nchar(get.ending(x, ending.type))>0, "_", ""), 
                                                                                   get.ending(x, ending.type))))
    }
  }
  study.data$UniqueLabel[grep("dilute", study.data$testSampleLabel)] = paste(study.data$UniqueLabel[grep("dilute", study.data$testSampleLabel)], "D", sep="_")
  if (ending.type == "all"){
    for( i in 1:nrow(study.data)){
      study.data$UniqueLabel[i] = ifelse(study.data$testSequencingPlatform[i]%in%Illumina.seq, paste0(study.data$UniqueLabel[i], "_I"), paste0(study.data$UniqueLabel[i], "_B"))
      study.data$UniqueLabel[i] = ifelse(length(grep("3300", study.data$testSampleLabel[i]))>0,
                                         paste0(study.data$UniqueLabel[i], "FF"), # fresh frozen sample
                                         ifelse(length(grep("5500", study.data$testSampleLabel[i]))>0,
                                                paste0(study.data$UniqueLabel[i], "DQ"), # diff quik samples
                                                ifelse(length(grep("6600", study.data$testSampleLabel[i]))>0, 
                                                       paste0(study.data$UniqueLabel[i], "CT"), # ctDNA samples
                                                       paste0(study.data$UniqueLabel[i], "FF")))) # the first of the fresh frozen samples before naming convention of 3300
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
  # Sample Label
  study.data$UniqueLabel = gsub(" ", "_", gsub("-", "_", study.data$UniqueLabel, fixed=TRUE))
  return(study.data)
}

# ***************
# function to filter methylation data for 450k & EPIC array imported from idat files
# ***************
filter.data = function(import.data, 
                       filter.type=c("filter","message"),
                       first.filter = TRUE,
                       array.type=array_type){
  if (first.filter) {
    meth.filter = champ.filter(beta=import.data$beta,
                               M=import.data$M,
                               pd=import.data$pd,
                               intensity=import.data$intensity,
                               Meth=import.data$Meth,
                               UnMeth=import.data$UnMeth,
                               detP=import.data$detP,
                               beadcount=import.data$beadcount,
                               arraytype=array_type,
                               filterNoCG=FALSE,
                               filterSNPs= FALSE,
                               filterMultiHit=FALSE,
                               filterXY=FALSE
                               )
  } else {
    meth.filter = champ.filter(beta=import.data$beta,
                               M=import.data$M,
                               pd=import.data$pd,
                               intensity=import.data$intensity,
                               Meth=import.data$Meth,
                               UnMeth=import.data$UnMeth,
                               detP=NULL,
                               beadcount=NULL,
                               arraytype=array_type,
                               autoimpute=FALSE,
                               filterDetP=FALSE,
                               filterBeads=FALSE,
                               filterSNPs= ifelse (array_type == "EPIC" , FALSE, TRUE)
                               )
  }
  if (!first.filter){
    if (array_type == "EPIC"){
      # *************
      # Filter EPIC SNPs using Zhou hg19 and McCartney paper
      # *************
      message("     Filtering SNPs from Zhou et.al. and McCartney et al. SNPS and multi-mappers.")
      before_num = nrow(meth.filter$beta)
      remove_probes = union(Zhou_SNP$IlmnID, union(McCartney_SNP$IlmnID, McCartney_multi$IlmnID))
      meth.filter = lapply(meth.filter, function(x) {if (length(colnames(x)) == length(import.data$pd$Sample_Name)){
        x[-which(rownames(x)%in%remove_probes),]} else {x}})
      removed_num = before_num - nrow(meth.filter$beta)
      message("    Removing ", removed_num, " probes from the analysis.")
    }
  }
  if (filter.type=="filter"){
    return(meth.filter)
  }
}

# ***************
# function to format the number reported in plots and tables
# ***************
format.number = function (number){
  num_char = ifelse(nchar(format(number, scientific=FALSE))>=7,
                    paste0(substr(format(number, scientific=FALSE), 1, nchar(format(number, scientific=FALSE))-6), ",",
                           substr(format(number, scientific=FALSE), nchar(format(number, scientific=FALSE))-5, nchar(format(number, scientific=FALSE))-3), ",",
                           substr(format(number, scientific=FALSE), nchar(format(number, scientific=FALSE))-2,
                                  nchar(format(number, scientific=FALSE)))),
                    ifelse(nchar(format(number, scientific=FALSE))>=4,
                           paste0(substr(format(number, scientific=FALSE), 1, nchar(format(number, scientific=FALSE))-3), ",",
                                  substr(format(number, scientific=FALSE), nchar(format(number, scientific=FALSE))-2,
                                         nchar(format(number, scientific=FALSE)))),
                           format(number, scientific=FALSE)))
  return(num_char)
}

# ***************
# function to extract the message in each step of the filtering process
# ***************
find.message = function(message.file, 
                        message.type, 
                        array.type,
                        norm.methods=normalisation_methods){
  if (length(grep("FunctionalNormalization", norm.methods))>0){
    switch(message.type,
           meth_pvalue = substr(message.file[grep("detection p-value above 0.01", message.file)], 
                                nchar(message.file[grep("detection p-value above 0.01", message.file)])/3,
                                nchar(message.file[grep("detection p-value above 0.01", message.file)])),
           meth_bead = substr(message.file[grep("beadcount <3 in at least 5% of samples", message.file)],
                              nchar(message.file[grep("beadcount <3 in at least 5% of samples", message.file)])-26,
                              nchar(message.file[grep("beadcount <3 in at least 5% of samples", message.file)])),
           meth_CpG = substr(message.file[grep("Filtering NoCG Start", message.file)+1], 22,
                             nchar(message.file[grep("Filtering NoCG Start", message.file)+1])),
           meth_XY = substr(message.file[grep("probes located on X,Y chromosome", message.file)],50,
                            nchar(message.file[grep("probes located on X,Y chromosome", message.file)])),
           meth_SNP = ifelse(array.type=="450K",
                             message.file[grep("Filtering SNPs Start", message.file)+3],
                             as.numeric(substr(message.file[grep("Filtering SNPs from Zhou", message.file)+1], 14, 19))+
                               as.numeric(substr(message.file[grep("Filtering SNPs Start", message.file)+3], 14, 19))),
           meth_multi = message.file[grep("Filtering MultiHit Start", message.file)+2],
           meth_samp = message.file[grep("samples remained for analysis.", message.file)]
    )
  } else {
    switch(message.type,
           meth_pvalue = substr(message.file[grep("detection p-value above 0.01", message.file)+1], 1,
                                nchar(message.file[grep("detection p-value above 0.01", message.file)+1])-1),
           meth_bead = message.file[grep("beadcount <3 in at least 5% of samples", message.file)+1],
           meth_CpG = substr(message.file[grep("Filtering NoCG Start", message.file)+1], 22,
                             nchar(message.file[grep("Filtering NoCG Start", message.file)+1])),
           meth_XY = substr(message.file[grep("probes located on X,Y chromosome", message.file)],50,
                            nchar(message.file[grep("probes located on X,Y chromosome", message.file)])),
           meth_SNP = ifelse(array.type=="450K",
                             message.file[grep("Filtering SNPs Start", message.file)+3],
                             message.file[grep("Filtering SNPs from Zhou", message.file)+1]),
           meth_multi = message.file[grep("Filtering MultiHit Start", message.file)+2],
           meth_samp = message.file[grep("samples remained for analysis.", message.file)]
    )
  }
}

# **********
# Set up sample_info with the sentrix id (combined slide and array) 
# **********
set.sample.info.levels = function(pd,
                                  samp.vars = sample_sheet_cat_all,
                                  variables = sample_categories, 
                                  var.levs = sample_levels){
  samp.info = pd
  if (!length(grep("sentrix", colnames(samp.info), ignore.case=TRUE))>0) {
    samp.info$SentrixID = vector(mode="character", length=nrow(samp.info))
    if (length(grep("slide", colnames(samp.info), ignore.case=TRUE))==1) {
      col.name.slide = grep("slide", colnames(samp.info), ignore.case=TRUE, value=TRUE)
      samp.info[,col_name_slide] = as.character(samp.info[,col.name.slide])
      } else  {
        col.name.slide = "Slide"
        samp.info[,col.name.slide] = vector(mode="character", length=nrow(samp.info))
        samp.info[,col.name.slide] = "No_Slide"
      }
    if (length(grep("array", colnames(samp.info), ignore.case=TRUE))==1) {
      col.name.array = grep("array", colnames(samp.info), ignore.case=TRUE, value=TRUE)
      samp.info[,col.name.array] = as.character(samp.info[,col.name.array])
      } else {
        col.name.array = "Array"
        samp.info[,col.name.array] = vector(mode="character", length=nrow(samp.info))
        samp.info[,col.name.array] = "No_Array"
      }
    samp.info$SentrixID = paste0(samp.info[,col.name.slide], "_", samp.info[,col.name.array])
    } else {
      sentrix.col = grep("sentrix", colnames(samp.info), ignore.case=TRUE)
      colnames(samp.info)[sentrix.col] = "SentrixID"
    }
  for (i in 1:length(variables)){
    var = variables[i]
    colnames(samp.info)[which(colnames(samp.info)==samp.vars[i])] = var
    if (length(grep("NA", var.levs[[var]]))>0){
      samp.info[,var][which(is.na(samp.info[,var]))] = "NA"
    }
    samp.info[,var] = droplevels(factor(samp.info[,var], levels=var.levs[[var]]))
    var.levs[[var]] = levels(samp.info[,var])
  }
  return(list(sample.info=samp.info, var.levels=var.levs))
}

# ***************
# function to format message and extract number in each step of the filtering process
# ***************
get.probe.info = function(message.file, 
                          message.type, 
                          array.type){
  if (is.list(message.file)){
    matches = 0
    for (i in 1:length(message.file)){
      probes.text = find.message(message.file[[i]], message.type, array.type)
      if (length(probes.text)>0){
        matches = matches + as.numeric(regmatches(probes.text, regexpr("[[:digit:]]+", probes.text)))
      } else {
        matches = matches + 0
      }
    }
  } else {
    probes.text = find.message(message.file, message.type, array.type)
    if (length(probes.text)>0){
      matches = regmatches(probes.text, regexpr("[[:digit:]]+", probes.text))
    } else {
      matches = 0
    }
  }
  if (message.type=="meth_samp"){
    filt.message = paste(format.number(as.numeric(matches)), "samples removed from the analysis.", sep=" ")
    } else {
      filt.message = paste(format.number(as.numeric(matches)), "probes removed from the analysis.", sep=" ")
    }
  return.data = list(as.numeric(matches), filt.message)
  names(return.data) = c("Number", "Message")
  return(return.data)
}

# ***************
# function for QC with changes to suit project (champ.QC -> QC)
# ***************
QC = function (raw=TRUE, 
               beta=meth_import$beta, 
               pheno=meth_import$pd$Sample_Group, 
               Group=NULL,
               col.cat = colour_cat_border$Sample_Group[order(names(colour_cat_border$Sample_Group))],
               mdsPlot=FALSE, 
               densityPlot=TRUE, 
               dendrogram=FALSE, 
               PDFplot=FALSE,
               Rplot=TRUE, 
               Feature.sel="None", 
               resultsDir=results_dir, 
               norm.method=1)
  {
  message("[===========================]")
  message("[<<<<< ChAMP.QC START >>>>>>]")
  message("-----------------------------")
  message("champ.QC Results will be saved in ", resultsDir)
  message("[QC plots will be proceed with ", dim(beta)[1],
          " probes and ", dim(beta)[2], " samples.]\\n")
  if (min(beta, na.rm = TRUE) == 0) {
    beta[beta == 0] <- 1e-06
    message("[", length(which(beta == 0)), " Zeros detect in your dataset, will be replaced with 0.000001]\\n")
    }
  if (ncol(beta) != length(pheno))
    stop("Dimension of DataSet Samples, pheno and name must be the same. Please check your input.")
  message("<< Prepare Data Over. >>")
  if (norm.method>1){
    norm.report = paste(normalisation_methods[1:norm.method], collapse="|")
  } else {
    norm.report = normalisation_methods[1]
  }
  if (densityPlot) {
    plot.title = ifelse(raw, 
                        paste0("Density plot of raw data (", format.number(nrow(beta)), " probes)"),
                        paste0("Density plot of ", norm.report, " normalised data (", format.number(nrow(beta)), " probes)"))
    if (PDFplot){
      if (!file.exists(resultsDir))
        dir.create(resultsDir)
      file.name = ifelse(raw, "raw_density.pdf", "normalised_densityPlot.pdf")
    } else {
      file.name = ""
    }
    if (Rplot)
      densityPlot(beta, sampGroups = pheno, main = plot.title, pal = col.cat, xlab = "Beta Value", cex.main=0.9, cex.lab=0.8)
    if (PDFplot) {
      pdf(paste(resultsDir, file.name, sep = "/"),
          width = 7, height = 5)
      densityPlot(beta, sampGroups = pheno, main = plot.title, pal = col.cat, xlab = "Beta Value", , cex.main=0.9, cex.lab=0.8)
      dev.off()
      }
    message("<< Plot densityPlot Done. >>")
  }
  if (mdsPlot) {
    plot.title = ifelse(raw, 
                        paste0("MDS plot of raw data \n(", format.number(numPos), " most variable positions)"), 
                        paste0("MDS plot of ", norm.report, " normalised data \n(", format.number(numPos), " most variable positions)"))
    if (PDFplot){
      if (!file.exists(resultsDir))
        dir.create(resultsDir)
      file.name = ifelse(raw, "raw_mdsPlot.pdf", "normalised_mdsPlot.pdf")
    } else {
      file.name = ""
    }
    if (Rplot)
      mdsPlot(beta, numPositions = numPos, sampGroups = pheno, pch=19, main = plot.title,
              colnames(beta))
    if (PDFplot) {
      pdf(paste(resultsDir, file.name, sep = "/"),
          width = 6, height = 6)
      mdsPlot(beta, numPositions = numPos, sampGroups = pheno, pch=19, main = plot.title)
      dev.off()
      }
    message("<< plot mdsPlot Done. >>\\n")
    }
  if (dendrogram) {
    if (Feature.sel == "None") {
      if (! is.null(Group)){
        colnames(beta) = Group
        }
      message("< Dendrogram Plot Feature Selection Method >: No Selection, directly use all CpGs to calculate distance matrix.")
      hc <- hclust(dist(t(beta)))
      }
    else if (Feature.sel == "SVD") {
      message("< Dendrogram Feature Selection Method >: Use top SVD CpGs to calculate distance matrix.")
      SVD <- svd(beta)
      rmt.o <- EstDimRMT(beta - rowMeans(beta))
      M <- SVD$v[, 1:rmt.o$dim]
      message(colnames(beta))
      rownames(M) <- colnames(beta)
      colnames(M) <- paste("Component", c(1:rmt.o$dim))
      hc <- hclust(dist(M))
      }
    dend <- as.dendrogram(hc)
    MyColor <- sample_colours[1:length(table(pheno))]
    names(MyColor) <- names(table(pheno))
    labels_colors(dend) <- MyColor[pheno[order.dendrogram(dend)]]
    dend <- dend %>% set("labels_cex", 0.8)
    dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 0.4) %>% set("leaves_col", MyColor[pheno[order.dendrogram(dend)]])
    plot.title = ifelse(raw, paste("All samples before normalization (", format.number(nrow(beta)), " probes)"),
                        paste("All samples after ", norm.report, " normalization (", format.number(nrow(beta)), " probes)"))
    if (PDFplot){
      if (!file.exists(resultsDir))
        dir.create(resultsDir)
      file.name = ifelse(raw, "raw_SampleCluster.pdf", "normalised_SampleCluster.pdf")
    } else {
      file.name = ""
    }
    if (Rplot) {
      plot(dend, center = TRUE, main = plot.title)
      legend("topright", fill = MyColor, legend = names(MyColor))
      }
    if (PDFplot) {
      pdf(paste(resultsDir, file.name, sep = "/"),
          width = floor(log(ncol(beta)) * 3), height = floor(log(ncol(beta)) *
                                                               2))
      plot(dend, center = TRUE, main = plot.title)
      legend("topright", fill = MyColor, legend = names(MyColor))
      dev.off()
      }
    message("<< Plot dendrogram Done. >>\\n")
    }
  message("[<<<<<< ChAMP.QC END >>>>>>>]")
  message("[===========================]")
  msg = ifelse(raw, "[You may want to process champ.norm() next.]\\n", "End of QC of normalised data.\\n")
  message(msg)
}

# ***************
# function for SVD with changes to suit project (champ.SVD -> SVD)
# ***************
SVD = function (beta=meth_norm, 
                rgSet=NULL, 
                pd=meth_import$pd, 
                RGEffect=FALSE, 
                PDFplot=FALSE, 
                Rplot=TRUE, 
                resultsDir="./CHAMP_SVDimages/", 
                batch=FALSE,
                type = "Methylation",
                x.text.size=0.8, 
                y.text.size=0.8) {
  message("[===========================]")
  message("[<<<<< ChAMP.SVD START >>>>>]")
  message("-----------------------------")
  GenPlot = function(thdens.o, estdens.o, evalues.v) {
    minx = min(min(thdens.o$lambda), min(evalues.v))
    maxx = max(max(thdens.o$lambda), max(evalues.v))
    miny = min(min(thdens.o$dens), min(estdens.o$y))
    maxy = max(max(thdens.o$dens), max(estdens.o$y))
    }
  EstDimRMTv2 = function(data.m) {
    M = data.m
    for (c in 1:ncol(M)) {
      M[, c] = (data.m[, c] - mean(data.m[, c]))/sqrt(var(data.m[, c]))
    }
    sigma2 = var(as.vector(M))
    Q = nrow(data.m)/ncol(data.m)
    thdens.o = thdens(Q, sigma2, ncol(data.m))
    C = 1/nrow(M) * t(M) %*% M
    eigen.o = eigen(C, symmetric = TRUE)
    estdens.o = density(eigen.o$values, from = min(eigen.o$values), to = max(eigen.o$values), cut = 0)
    GenPlot(thdens.o, estdens.o, eigen.o$values)
    intdim = length(which(eigen.o$values > thdens.o$max))
    return(list(cor = C, dim = intdim, estdens = estdens.o, thdens = thdens.o))
  }
  
  thdens = function(Q, sigma2, ns) {
    lambdaMAX = sigma2 * (1 + 1/Q + 2 * sqrt(1/Q))
    lambdaMIN = sigma2 * (1 + 1/Q - 2 * sqrt(1/Q))
    delta = lambdaMAX - lambdaMIN
    roundN = 3
    step = round(delta/ns, roundN)
    while (step == 0) {
      roundN = roundN + 1
      step = round(delta/ns, roundN)
      }
    lambda.v = seq(lambdaMIN, lambdaMAX, by = step)
    dens.v = vector()
    ii = 1
    for (i in lambda.v) {
      dens.v[ii] = (Q/(2 * pi * sigma2)) * sqrt((lambdaMAX - i) * (i - lambdaMIN))/i
      ii = ii + 1
      }
    return(list(min = lambdaMIN, max = lambdaMAX, step = step, lambda = lambda.v, dens = dens.v))
  }
  
  message("champ.SVD Results will be saved in ", resultsDir, " .\\n")
  if (length(which(is.na(beta))) > 0){
    message(length(which(is.na(beta))), " NA are detected in your beta Data Set, which may cause failure of or incorrect SVD analysis. Impute NA with champ.impute() function first.")
  }
  message("[SVD analysis will be proceed with ", dim(beta)[1], " probes and ", dim(beta)[2], " samples.]\\n")
  message("\\n[ champ.SVD() will only check the dimensions between data and pd.\\n")
  message("Sample_Names must be correctly matched between data and pd file), thus please ensure pd file sample names is in accordance with data sets (beta and rgSet).\\n")
  if (is.null(pd) | class(pd) == "list") {
    stop("pd parameter in Data Frame or Matrix is necessary And must contain at least two factors. If your pd is a list, please change its Format.")
  }
  if (class(pd) == "matrix"){
    pd = as.data.frame(pd)
  }
  PhenoTypes.lv_tmp = pd[, !colnames(pd) %in% c("Sample_Name", "Project", "filenames", "Basename") & apply(pd, 2, function(x) length(unique(x))) != 1]
  PhenoTypes.lv = PhenoTypes.lv_tmp
  for (i in ncol(PhenoTypes.lv_tmp)) {
    if (class(PhenoTypes.lv_tmp[, i]) != "numeric") 
      PhenoTypes.lv[, i] = as.factor(as.numeric(as.factor(PhenoTypes.lv_tmp[, i])))
  }
  
  if (!is.null(rownames(pd))) 
    rownames(PhenoTypes.lv) = rownames(pd)
  if (ncol(PhenoTypes.lv) >= 2) {
    message("<< Following Factors in your pd(sample_sheet.csv) will be analysised: >>")
    sapply(colnames(PhenoTypes.lv_tmp), function(x) message("<", x, ">(", class(PhenoTypes.lv[[x]]), "):", paste(unique(PhenoTypes.lv_tmp[, x]), collapse = ", ")))
    message("champ.SVD automatically selects ALL factors containing at least two different values from your pd(sample_sheet.csv)\\n")
    message("Please remove manually from your pd any variable not part of the analysis, then retry champ.SVD().")
  } else {
    stop("You don't a have a factor with at least two values in analysis. If factors contain only one value, no variation with be found...")
  }
  
  if (ncol(pd) > ncol(PhenoTypes.lv)) {
    message("\\n<< Following Factors in your pd(sample_sheet.csv) will not be in analysis: >>")
    sapply(setdiff(colnames(pd), colnames(PhenoTypes.lv)), function(x) message("<", x, ">"))
    message("[Factors are ignored because they only indicate Name or Project, or they contain ONLY ONE value across all Samples.]")
  }
  
  if (RGEffect == TRUE & is.null(rgSet)) {
    message("If you want to check Effect of Control Probes, you MUST provide rgSet parameter. Now champ.SVD can only analysis factors in pd.")
  }
  if (!is.null(rgSet) & RGEffect) {
    if (rgSet@annotation[1] == "IlluminaHumanMethylation450k"){
      data(ControlProbes450K)
    } else {
      data(ControlProbesEPIC)
    }
    dataC2.m <- as.data.frame(log2(apply(ControlProbes, 1, function(x) if (x[3] == "Grn") getGreen(rgSet)[x[2],] else getRed(rgSet)[x[2], ])))
    PhenoTypes.lv = cbind(PhenoTypes.lv, dataC2.m)
    message("\\n<< Following rgSet information have been added to PhenoTypes.lv. >>")
    sapply(colnames(dataC2.m), function(x) message("<", x, ">")) 
    }
  if (nrow(PhenoTypes.lv) == ncol(beta)){
    message("\\n<< PhenoTypes.lv generated successfully. >>")
  } else {
    stop("Dimension of your pd file (and rgSet information) is not equal to your beta matrix.")
  }
  tmp.m = beta - rowMeans(beta)
  rmt.o = EstDimRMTv2(tmp.m)
  svd.o = svd(tmp.m)
  if (rmt.o$dim > 20){
    topPCA = 20
  } else {
    topPCA = rmt.o$dim
  }
  svdPV.m = matrix(nrow = topPCA, ncol = ncol(PhenoTypes.lv))
  colnames(svdPV.m) = colnames(PhenoTypes.lv)
  for (c in 1:topPCA) {
    for (f in 1:ncol(PhenoTypes.lv)) {
      if (class(PhenoTypes.lv[, f]) != "numeric") {
        svdPV.m[c, f] = kruskal.test(svd.o$v[, c] ~ as.factor(PhenoTypes.lv[[f]]))$p.value
      } else {
        svdPV.m[c, f] = summary(lm(svd.o$v[, c] ~ PhenoTypes.lv[[f]]))$coeff[2, 4]
      }
    }
  }
  
  message("<< Calculate SVD matrix successfully. >>")
  myPalette = c(brewer.pal(9, "YlOrRd")[9], brewer.pal(9, "YlOrRd")[8], "#FD4321", brewer.pal(9, "YlOrRd")[5],
                brewer.pal(9, "YlOrRd")[4], brewer.pal(9, "YlOrRd")[3], brewer.pal(9, "YlOrRd")[2], "#FEFEF1")
  breaks.v = c(-200, -20, -10, -7, -5, -3, -1, log10(0.05), 0)
  if (Rplot) {
    par(mar = c(5, 10, 4, 3))
    batch_msg = ifelse(batch, " with batch correction.", " without batch correction.")
    normalisation = ifelse(type=="Methylation", paste(normalisation_methods, collapse="|"), "RNASeq TMM")
    image(x = 1:nrow(svdPV.m), y = 1:ncol(svdPV.m), z = log10(svdPV.m), 
          col = myPalette, breaks = breaks.v, xlab = "", ylab = "", 
          axes = FALSE, 
          main = paste0("Singular Value Decomposition Analysis (SVD)\n", 
                       normalisation, " normalised data", batch_msg))
    axis(1, at = 1:nrow(svdPV.m), labels = paste0("PC-", 1:nrow(svdPV.m)), las = 2, cex.axis=x.text.size)
    suppressWarnings(axis(2, at = 1:ncol(svdPV.m), labels = colnames(svdPV.m), las = 2, cex.axis=y.text.size))
    legend("topright", #x = -(topPCA/2.5), y = 3, 
           legend = c(expression("p < 1x" ~ 10^{-20}), expression("p < 1x" ~ 10^{-10}), expression("p < 1x" ~ 10^{-7}),
                      expression("p < 1x" ~ 10^{-5}), expression("p < 1x" ~ 10^{-3}), "p < 0.01", "p < 0.05", "p > 0.05"), 
           cex=0.9, fill = myPalette, col="black", bty="n")
  }
  
  if (PDFplot) {
    if (!file.exists(resultsDir))
        dir.create(resultsDir)
    pdf(paste0(resultsDir, "SVDsummary.pdf"), width = 8, 
        height = 8)
    par(mar = c(5, 10, 2, 1), xpd = TRUE)
    image(x = 1:nrow(svdPV.m), y = 1:ncol(svdPV.m), z = log10(svdPV.m), 
          col = myPalette, breaks = breaks.v, xlab = "", ylab = "", 
          axes = FALSE, main = paste0("Singular Value Decomposition Analysis (SVD)\n", paste(normalisation_methods, collapse="|"), " normalisation", batch_msg))
    axis(1, at = 1:nrow(svdPV.m), labels = paste0("PC-", 1:nrow(svdPV.m)), las = 2)
    suppressWarnings(axis(2, at = 1:ncol(svdPV.m), labels = colnames(svdPV.m), las = 2))
    legend(x = -(topPCA/2.5), y = 3, 
           legend = c(expression("p < 1x" ~ 10^{-10}), expression("p < 1x" ~ 10^{-5}), "p < 0.01", "p < 0.05", "p > 0.05"), 
           fill = myPalette,  pch=22, cex=0.5, pt.cex=0.8, bty="n", par("usr")[2], par("usr")[4], xpd = NA)
    dev.off()
  }
  
  message("<< Plot SVD matrix successfully. >>")
  message("[<<<<<< ChAMP.SVD END >>>>>>]")
  message("[===========================]")
  message("[If the batch effect is not significant, you may want to process champ.DMP() or champ.DMR() or champ.BlockFinder() next]")
  message("[otherwise, you may want to run champ.runCombat() to eliminat batch effect, then rerun champ.SVD() to check corrected result.]")
  return(svdPV.m)
}

# ***************
# function for runCombat with changes to suit project (champ.runCombat -> runCombat)
# ***************
runCombat = function (beta=meth_norm, 
                      pd=meth_import$pd,
                      variablename="Sample_Group", 
                      batchname=c("Slide"), 
                      logitTrans = TRUE) 
  {
  message("[===========================]")
  message("[<< CHAMP.RUNCOMBAT START >>]")
  message("-----------------------------")
  message("<< Preparing files for ComBat >>")
  if (length(which(is.na(beta))) > 0) 
    message(length(which(is.na(beta))), " NA are detected in your beta Data Set, which may cause fail or uncorrect of runCombat analysis. You may want to impute NA with champ.impute() function first.")
  message("[Combat correction will be proceed with ", dim(beta)[1], 
          " probes and ", dim(beta)[2], " samples.]\\n")
  if (is.null(pd) | class(pd) == "list") 
    stop("pd parameter in Data Frame or Matrix is necessary And must contain at least tow factors. If your pd is a list, please change its Format.")
  if (class(pd) == "matrix") 
    pd <- as.data.frame(pd)
  if (is.null(variablename) | !variablename %in% colnames(pd)) 
    stop("variablename parameter MUST contains variable in pd file.")
  valid.idx <- which(!colnames(pd) == variablename & apply(pd, 2, function(x) length(unique(x))) != 1 & apply(pd, 2, function(x) all(table(x) >= 2)))
  if (length(valid.idx) == 0) 
    stop("There is no valid factor can be corrected. Factor can be corrected must contian at least two phenotypes, each of them must contain at least two samples. Also batch factors can be variable factor.
         Please check if your covariates fulfill these requirement.")
  PhenoTypes.lv_tmp <- as.data.frame(pd[, valid.idx])
  colnames(PhenoTypes.lv_tmp) <- colnames(pd)[valid.idx]
  PhenoTypes.lv <- as.data.frame(apply(PhenoTypes.lv_tmp, 2, function(x) if (class(x) != "numeric") as.factor(as.numeric(as.factor(x)))))
  if (!is.null(rownames(pd))) 
    rownames(PhenoTypes.lv) <- rownames(pd)
  if (ncol(PhenoTypes.lv) >= 1) {
    message("<< Following Factors in your pd(sample_sheet.csv) could be applied to Combat: >>")
    sapply(colnames(PhenoTypes.lv_tmp), function(x) message("<", 
                                                            x, ">(", class(PhenoTypes.lv[[x]]), "):", paste(unique(PhenoTypes.lv_tmp[, x]), collapse = ", ")))
    message("[champ.runCombat have automatically select ALL factors contain at least two different values from your pd(sample_sheet.csv).]")
    }
  else {
    stop("You don't have even one factor with at least two value to be analysis. Maybe your factors contains only one value, no variation at all...")
    }
  if (ncol(pd) > ncol(PhenoTypes.lv)) {
    message("\\n<< Following Factors in your pd(sample_sheet.csv) can not be corrected: >>")
    sapply(setdiff(colnames(pd), colnames(PhenoTypes.lv)), 
           function(x) message("<", x, ">"))
    message("[Factors are ignored because they are conflict with variablename, or they contain ONLY ONE value across all Samples, or some phenotype contains less than 2 Samples.]")
    }
  if (all(batchname %in% colnames(PhenoTypes.lv))) {
    message("As your assigned in batchname parameter: ", 
            paste(batchname, collapse = ","), " will be corrected by Combat function.")
    }
  else {
    stop(setdiff(batchname, colnames(PhenoTypes.lv)), " factors is not valid to run Combat, please recheck your dataset.")
    }
  for (i in 1:length(batchname)) {
    message("\\n<< Checking confounded status between ", batchname[i], 
            " and ", variablename, " >>")
    if (i + 1 <= length(batchname)) 
      tmp_formdf <- as.formula(paste(" ~ ", paste(c(variablename, 
                                                    batchname[(i + 1):length(batchname)]), collapse = " + "), 
                                     sep = ""))
    else tmp_formdf <- as.formula(paste(" ~", variablename))
    message("--------------------------")
    message("Model for Correction is: ")
    print(tmp_formdf)
    tmp_mod <- model.matrix(tmp_formdf, data = pd)
    tmp_batchmod <- model.matrix(~-1 + pd[, batchname[i]])
    n.batch <- length(unique(pd[, batchname[i]]))
    tmp_design <- cbind(tmp_batchmod, tmp_mod)
    check <- apply(tmp_design, 2, function(x) all(x == 1))
    tmp_design <- as.matrix(tmp_design[, !check])
    message("Combat can adjust for ", ncol(tmp_design) - ncol(tmp_batchmod), 
            " covariate(s) or covariate level(s)\\n")
    if (qr(tmp_design)$rank < ncol(tmp_design)) {
      message("\\nSorry, Your covariate is confounded with batch, which means, some of your covariates or can be totally represented by batches. In other word, your Model generated from ", 
              variablename, " and ", batchname[i], " for correction is not a full rank matrix, which is not in accord of the condition of Combat.\\n")
      stop("Sorry, your covariate ", variablename, " is confounded with batch ", 
           batchname[i])
      }
    message("--------------------------")
    }
  message("<< Rank Check Complete, you data is good to proceed. >> ^_^")
  if (min(beta) <= 0) {
    message("Zeros in your dataset have been replaced with smallest positive value.")
    beta[beta <= 0] <- min(beta[beta > 0])
    }
  beta_2 <- beta
  innercombat <- function(beta, batch, formdf) {
    if (logitTrans) 
      beta <- logit2(beta)
    print(formdf)
    mod <- model.matrix(formdf, data = pd)
    message("Generate mod success. Started to run ComBat, which is quite slow...")
    library("sva", lib.loc="/software/R/R-3.4.1/lib64/R/library")
    combat <- ComBat(dat = beta, batch = batch, mod = mod, 
                     par.prior = TRUE)
    if (logitTrans) 
      combat = ilogit2(combat)
    return(combat)
    }
  for (i in 1:length(batchname)) {
    message("\\n<< Start Correcting ", batchname[i], " >>")
    if (i + 1 <= length(batchname)) 
      formdf <- as.formula(paste(" ~ ", paste(c(variablename, 
                                                batchname[(i + 1):length(batchname)]), collapse = " + "), 
                                 sep = ""))
    else formdf <- as.formula(paste(" ~", variablename))
    beta <- innercombat(beta, pd[, batchname[i]], formdf)
    }
  if (is.null(beta)) {
    message("Sorry, champ.runCombat failed. Your old dataset will be returned.")
    return(beta_2)
    }
  else {
    message("champ.runCombat success. Corrected dataset will be returned.")
    return(beta)
    }
  message("[<<< CHAMP.RUNCOMBAT END >>>]")
  message("[===========================]")
  message("[You may want to process champ.SVD() next.]\\n")
}

# ***************
# function to normalise data using one or combination of "QN", "BMIQ", "PBC"
# ***************
normalise.data = function(dat, 
                          rgSet,
                          norm.method, 
                          array.type,
                          seed){
  set.seed(seed)
  switch(norm.method,
         QN = normalize.quantiles.robust(dat, remove.extreme="none"), 
         BMIQ = champ.norm(beta=dat, 
                           method=norm.method, 
                           resultsDir=paste0(results_dir, norm.method, "_normalisation/"), 
                           plotBMIQ=FALSE, 
                           arraytype=array.type,
                           cores=1),
         PBC = champ.norm(beta=dat, 
                          method=norm.method, 
                          resultsDir=paste0(results_dir, norm.method, "_normalisation/"), 
                          plotBMIQ=FALSE, 
                          arraytype=array.type,
                          cores=1),
         FunctionalNormalization = champ.norm(beta=dat,
                                              rgSet=rgSet,
                                              method=norm.method,
                                              resultsDir=paste0(results_dir, norm.method, "_normalisation/"), 
                                              plotBMIQ=FALSE, 
                                              arraytype=array.type,
                                              cores=1))
}

# ***************
# Function to look at the fit of dendogram cluster
# ***************
get.dend.fit = function(logit.mat,
                        clust.meth = "ward.D",
                        dist.meth = "euclidean",
                        k.values){
  # bootstrap clustering methods
  eval.df = data.frame(Clusters = c(1:max(k.values), "Avg"),
                       Best_fit = vector(mode="character", length=max(k.values)+1))
  for (k.val in k.values){
    message(paste("Running bootstrap for", clust.meth, "with", k.val, "clusters"))
    fit = clusterboot(data = t(logit.mat), B=100, multipleboot = TRUE, clustermethod = hclustCBI, method=clust.meth, scaling=FALSE, k = k.val, count=FALSE, seed=100) 
    eval.df = cbind(eval.df, c(fit$bootmean, rep(0,max(k.values)-k.val), mean(fit$bootmean)), c(fit$bootbrd/100, rep(0,max(k.values)-k.val), mean(fit$bootbrd/100)), c(unlist(lapply(fit$result$clusterlist, function(x) length(which(x==TRUE)))), rep(0,max(k.values)-k.val+1)))
    colnames(eval.df)[c((ncol(eval.df)-2):ncol(eval.df))] = c(paste("AvgJaccard",k.val, sep="_"), paste("Instability",k.val, sep="_"), paste("Num_in_Clust", k.val, sep="_"))
  }
  eval.df$Best_fit = apply(eval.df[, grep("AvgJaccard", colnames(eval.df))], 1, function(x) gsub("AvgJaccard_", "", colnames(eval.df[,grep("AvgJaccard", colnames(eval.df))])[which(x==max(x))[1]]))
  
  title = textGrob("Summary of Fit for Clustering Methods.", gp=gpar(fontsize=12))
  padding = unit(5, "mm")
  colnames(eval.df) = gsub("AvgJaccard_", "Fit ", gsub("Instability ", "Instblty ", colnames(eval.df)))
  eval.table = gtable_add_rows(tableGrob(eval.df,
                                         theme=tt1, row=NULL),
                               heights = grobHeight(title) + padding,
                               pos = 0)
  eval.table = gtable_add_grob(
    eval.table,
    title,
    1, 1, 1, ncol(eval.table)
  )
  return(eval.table)
}

# ***************
# Function to produce the sample dendograms on the dichotomised data
# ***************
# cluster the tumour samples using binary distance metric and Ward's method for linkage
sample.dendogram = function (dat,
                             samp.type = c("tumour", "all"),
                             probe.type = "",
                             clust.col,
                             k.num = 4,
                             col.lab.cex = 0.4,
                             main.cex = 1,
                             branch.wdth = 2,
                             dist.method = "binary",
                             clust.method = "ward.D"
){
  dend.meth = t(dat) %>% dist(method=dist.method) %>% hclust(method=clust.method) %>% as.dendrogram
  if(samp.type=="tumour"){
    dend.colours = clust.col[2:(k.num+1)]
  } else {
    dend.colours = clust.col[1:k.num]
  }
  dend.meth %>%
    color_branches(k=k.num, col=dend.colours) %>%
    set("labels_cex", col.lab.cex) %>%
    set("branches_lwd", branch.wdth) %>%
    plot(main=paste0(study_name, " - ", paste0(ifelse(samp.type=="tumour", "Tumour ", "All")), 
                    " Samples \n", format.number(nrow(dat)), " ", probe.type, 
                    " probes in ", k.num, " Clusters."), 
         cex.main=main.cex)
  dend.heights = (unique(dend.meth %>% get_nodes_attr("height")))[order(unique(dend.meth %>% get_nodes_attr("height")), decreasing=TRUE)]
  v.height = (dend.heights[k.num-1] + dend.heights[k.num])/2
  abline(h=v.height, lty=2, lwd=2, col="red")
}

# ***************
# Function to extract the cluster and order of a dendogram on the dichotomised data
# ***************
get.leaf.order = function(dat,
                          dist.method="binary",
                          hclust.method="ward.D",
                          k.num=4){
  # get the subgroup classifications for the samples
  dend.meth = t(dat) %>% dist(method=dist.method) %>% hclust(method=hclust.method) %>% as.dendrogram
  cut.meth = t(dat) %>% dist(method=dist.method) %>% hclust(method=hclust.method)%>% cutree(k=k.num)
  cut.meth = cut.meth[labels(dend.meth)]
  # Labels subgroups in dendogram labels order
  leaf.clust = vector(mode="character", length=length(cut.meth))
  if (length(unique(cut.meth[grep("_T", names(cut.meth))]))>3){
    tumour.dend = t(dat[,grep("_T", colnames(dat))]) %>% dist(method=dist.method) %>% hclust(method=hclust.method)
    tumour.cut = t(dat[,grep("_T", colnames(dat))]) %>% 
      dist(method=dist.method) %>% 
      hclust(method=hclust.method) %>% cutree(k=k.num)
    tumour.cut = tumour.cut[labels(tumour.dend)]
    tumour.clusters = unique(tumour.cut)
    cluster.num = 1
    for (cluster in tumour.clusters){
      tumour.cut[which(tumour.cut==cluster)] = paste0("Group", cluster.num)
      cluster.num = cluster.num + 1
    }
    leaf.clust = cut.meth[grep("_T", names(cut.meth), invert=TRUE, fixed=TRUE, value=TRUE)] 
    leaf.clust = c(leaf.clust, tumour.cut)
    names(leaf.clust) = c(grep("_T", names(cut.meth), invert=TRUE, value=TRUE, fixed=TRUE), names(tumour.cut))
    leaf.clust[grep("_T", names(leaf.clust), invert=TRUE, fixed=TRUE)] = "Norm"
  } else {
    names(leaf.clust) = names(cut.meth)
    # Rename Clusters to reflect normal and tumour subtypes
    leaf.clust[grep("_T", names(cut.meth), invert=TRUE, fixed=TRUE)] = "Norm"
    # Extract unique tumour subgroups 
    tumour.clusters = unique(cut.meth[grep("_T", names(cut.meth), fixed=TRUE)])
    cluster.num = 1
    for (cluster in tumour.clusters){
      leaf.clust[which(cut.meth==cluster)] = paste0("Group", cluster.num)
      cluster.num = cluster.num + 1
    }
  }
  return(leaf.clust)
}

# ***************
# Create function to create a probe order based on copy number clustering (i.e. by chromosome arms)
# ***************
get.probe.order = function(dat,
                           probe.info=Probe_info,
                           cent.info=centromere_info,
                           probe.chr.levels){
  probe.names = NULL
  for (i in 1:22){
    new.probes = dplyr::filter(probe.info, CHR==as.character(i), IlmnID%in%rownames(dat))
    new.probes %<>% arrange(MAPINFO) %>% mutate(CHR_sub = ifelse(MAPINFO<cent.info$Centromere_start[i], 
                                                                 paste0(i, "p"),
                                                                 ifelse(MAPINFO>cent.info$Centromere_end[i], paste0(i, "q"), "0"))) %>% dplyr::select(IlmnID, CHR, MAPINFO, CHR_sub)
    probe.names = rbind(probe.names, new.probes)
  }
 probe.names$CHR_sub = as.character(probe.names$CHR_sub)
  start=TRUE
  for (i in 1:length(probe.names$CHR_sub)){
    if (start){
      new_chr = probe.names$CHR_sub[i]
      probe.names$CHR_sub[i] = paste("Chr:", new_chr, sep="")
      start = FALSE
    } else if (probe.names$CHR_sub[i] == new_chr){
      probe.names$CHR_sub[i] = " "
    } else {
      new_chr = probe.names$CHR_sub[i]
      probe.names$CHR_sub[i] = paste("Chr:", new_chr, sep="")
    }
  }
  return(probe.names)
}

# ***************
# Create function to produce the heatmaps from sample dendograms on the dichotomised and beta values
# ***************
sample.dend.heatmap = function(dat,
                               meth.norm = meth_norm,
                               samp.inf = sample_info,
                               samp.type = c("tumour", "normal", "all"),
                               samp.clust.col = samp_clust_col,
                               clust.rows=TRUE,
                               probe.levels = NULL,
                               Annotations=TRUE,
                               row.annotation = FALSE,
                               probe.type = "",
                               k.num = 4,
                               hmap.txt = hmap_txt_size,
                               hmap.col = hmap_col_size,
                               samp.cat = sample_categories,
                               col.cat = colour_cat,
                               clust.row.num = 2,
                               cutoff = cutoff,
                               silent=FALSE,
                               samp.clust=NULL,
                               dist.method = "euclidean",
                               hclust.method = "ward.D2",
                               new.meth.subroup = "Methylation_Subgroup"
                               ){
  if (length(grep("_T", colnames(meth.norm)))==0){
    colnames(meth.norm)[grep("benign", samp.inf$Sample_Group, ignore.case = TRUE, invert = TRUE)] = paste0(colnames(meth.norm)[grep("benign", samp.inf$Sample_Group, ignore.case = TRUE, invert = TRUE)], "_T")
  }
  # get the subgroup classifications for the samples
  if (is.null(samp.clust)){
    samp.clust = get.leaf.order (dat,
                                 dist.method=dist.method,
                                 hclust.method=hclust.method,
                                 k.num=k.num)
    if (length(grep("_B", names(samp.clust)))>0){
      # Rename Clusters to reflect normal and tumour subtypes
      samp.clust[grep("_B", names(samp.clust))] = "Bgn"
      # Extract unique tumour subgroups 
      tumour.clusters = unique(samp.clust[grep("_T", names(samp.clust))])
      cluster.num = 1
      for (cluster in tumour.clusters){
        samp.clust[which(samp.clust==cluster)] = paste0("Group", cluster.num)
        cluster.num = cluster.num + 1
      }
      samp.clust = samp.clust[order(samp.clust, decreasing = FALSE)]
    }
  }
  if (!clust.rows){
    probe.order = get.probe.order(dat,
                                  probe.info=Probe_info,
                                  cent.info=centromere_info,
                                  probe.chr.levels=probe.levels)
    probe.names = as.character(probe.order$IlmnID)
    dat.mat = dat[probe.names, names(samp.clust)]
    meth.norm.mat = meth.norm[probe.names,names(samp.clust)]
    rownames(dat.mat) = probe.order$CHR_sub
    rownames(meth.norm.mat) = probe.order$CHR_sub
  } else {
    dat.mat = dat[,names(samp.clust)]
    meth.norm.mat = meth.norm[which(rownames(meth.norm)%in%rownames(dat.mat)),
                              which(colnames(meth.norm)%in%colnames(dat.mat))]
    meth.norm.mat = meth.norm.mat[,names(samp.clust)]
  }
  norm.clust = samp.clust[grep("_B", names(samp.clust))]
  # create a dataframe of the sample categories to annotate at the top of the heatmap
  if (samp.type=="tumour"){
    if (is.null(samp.inf)){
      new.samp.inf = NULL
    } else{
      new.samp.inf = samp.inf[grep("_T", samp.inf$Sample_Name), ]
    }
    samp.clust = samp.clust[grep("_T", names(samp.clust))]
    meth.norm.mat = meth.norm.mat[ ,grep("_T", colnames(meth.norm.mat))]
  } else {
    new.samp.inf = samp.inf
  }
  if (!is.null(new.samp.inf)){
    rownames(new.samp.inf) = new.samp.inf$Sample_Name
    new.samp.inf = new.samp.inf[names(samp.clust),]
  }
  if (Annotations) {
    if (!is.null(new.samp.inf)){
      annotations = data.frame(New_subgroup=factor(samp.clust, levels=c(unique(samp.clust))),
                               new.samp.inf[, samp.cat]) 
      colnames(annotations) = c(new.meth.subroup, samp.cat)
      anno.cols = list()
      if (samp.type=="tumour"){
        anno.cols[[new.meth.subroup]] = samp.clust.col[2:(length(unique(samp.clust))+1)]
        names(anno.cols[[new.meth.subroup]]) = as.character(unique(samp.clust[grep("_T", names(samp.clust))]))
      } else {
        anno.cols[[new.meth.subroup]] = samp.clust.col[1:(length(unique(samp.clust)))]
        names(anno.cols[[1]]) = as.character(unique(samp.clust))
      }
      for (cat in names(samp.cat)){
        samp.levels = col.cat[[cat]]
        anno.cols[[cat]] = c(samp.levels)
      }
      names(anno.cols) = c(new.meth.subroup, samp.cat)
      rownames(annotations) = colnames(dat.mat)
    } 
  } else {
    annotations = data.frame(New_subgroup=factor(samp.clust, levels=c(unique(samp.clust))))
    colnames(annotations) = new.meth.subroup
    anno.cols = list()
    if (samp.type=="tumour"){
      anno.cols[[1]] = samp.clust.col[2:(length(unique(samp.clust))+1)]
      names(anno.cols[[1]]) = as.character(unique(samp.clust[grep("_T", names(samp.clust))]))
    } else {
      anno.cols[[1]] = samp.clust.col[1:(length(unique(samp.clust)))]
      names(anno.cols[[1]]) = as.character(unique(samp.clust))
    }
    names(anno.cols) = c(new.meth.subroup)
  }
  probes.to.report = format.number(nrow(dat.mat))
  # calculate normal tissue row annotation if required
  if (row.annotation){
    row.annotations = as.data.frame(apply(meth.norm[rownames(meth.norm.mat), names(norm.clust)], 2, function(x) as.character(round(x,2)*100)))
    rownames(row.annotations) = rownames(meth.norm.mat)
    row.cols = c(colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
    names(row.cols) = as.character(1:100)
    for (col in colnames(row.annotations)){
      row.annotations[,col] = factor(row.annotations[,col], levels = as.character(1:100))
      anno.cols[[col]] = row.cols
    }
  } else {
    row.annotations = NULL
  }
  # create the heatmap
  main.title = paste0(study_name, " - ", probes.to.report, " ", probe.type, " probes with ", dev_type, " >=", cutoff, ".\n")
  samp.clust = samp.clust[order(samp.clust, decreasing = FALSE)]
  if ("Bgn"%in%samp.clust){
    col.clust = length(grep("Bgn", samp_clust))
  } else {
    col.clust = c()
  }
  for (i in 1:(k.num-1)){
    col.clust = c(col.clust, length(grep(gsub(" ", "", paste("Bgn", paste("Group", 1:i, collapse="|"), sep="|")), samp.clust)))
  }
  row.clust = c()
  if (!clust.rows){
    row.pos = 0
    for (chr in unique(probe.order$CHR)[1:(length(unique(probe.order$CHR))-1)]){
      row.pos = row.pos + length(which(probe.order$CHR==chr))
      row.clust = c(row.clust, row.pos)
    }
  } else {
    row.clust = NULL
  }
  beta.hm = pheatmap::pheatmap(meth.norm.mat,   
                               annotation_col= annotations,
                               annotation_row = row.annotations,
                               scale="none",
                               clustering_method="ward.D2",
                               cluster_rows=clust.rows,
                               annotation_colors=anno.cols,
                               annotation_legend=ifelse(is.null(annotations), FALSE, TRUE),
                               annotation_names_row=ifelse(is.null(annotations), FALSE, TRUE),
                               show_rownames=ifelse(clust.rows, FALSE, TRUE),
                               fontsize=hmap.txt,
                               # cutree_cols=ifelse(samp.type=="tumour", k.num-1, k.num),
                               cutree_rows=clust.row.num*2,
                               fontsize_col=hmap.col,
                               cluster_cols=FALSE,
                               gaps_col = col.clust,
                               gaps_row = row.clust,
                               main=ifelse(!Annotations, "", paste0(main.title, " Supervised by Sample & Clustered by Probe beta values.")),
                               silent=silent
                               )
  # if using row annotations to add sample viewing change formating
  # LoH_dend_hmaps$beta$gtable$layout$name
  # [1] "main"                 "row_tree"             "matrix"               "col_names"            "col_annotation"       "col_annotation_names" "row_annotation"       "row_annotation_names"
  # [9] "annotation_legend"    "legend"   
  if (row.annotation){
    # [[8]] "row_annotation_names"
    beta.hm$gtable$grobs[[8]]$gp$fontsize = hmap.col # [[8]] "row_annotation_names"
    beta.hm$gtable$grobs[[8]]$x$arg2 = 0.8*beta.hm$gtable$grobs[[8]]$x$arg2
    #[[7]] "row_annotation"
    beta.hm$gtable$grobs[[7]]$x$arg2 = 0.8*beta.hm$gtable$grobs[[7]]$x$arg2 
    beta.hm$gtable$grobs[[7]]$width = unit(5, "bigpts")
    # change the overall width of the row annotations
    beta.hm$gtable$widths[[2]] = 0.85*beta.hm$gtable$widths[[2]]
    # remove sample annotation from legend [[9]] "legend"
    beta.hm$gtable$grobs[[9]]$children = beta.hm$gtable$grobs[[9]]$children[1:(ncol(annotations)*3)]
    beta.hm$gtable$grobs[[9]]$childrenOrder = beta.hm$gtable$grobs[[9]]$childrenOrder[1:(ncol(annotations)*3)]
    for (lab in names(beta.hm$gtable$grobs[[9]]$children)){ 
      if (!is.null(beta.hm$gtable$grobs[[9]]$children[[lab]])){
        beta.hm$gtable$grobs[[9]]$children[[lab]]$gp$fontsize=8
      }
    }
  }
  sample.heatmap = list(beta.hm)
  names(sample.heatmap) = c("beta")
  return(sample.heatmap)
}

# ***************
# A function to add the new sample subgroups to the sample_info
# ***************
add.samp.cat = function(cat.name=NULL,
                        cat.values=NULL,
                        cat.levels=NULL,
                        col=NULL,
                        col.bord=NULL,
                        samp.inf=sample_info,
                        samp.cat=sample_categories,
                        samp.type="",
                        samp.lev=sample_levels,
                        cat.cols=colour_cat,
                        cat.col.borders=colour_cat_border,
                        samp.clust = NULL
                        ){
  
  if (is.null(cat.name)||is.null(cat.values)||is.null(col)||is.null(col.bord)){
    message("<< Add sample category >>")
    message("This function must be supplied with the category name, vector of values, level colours and level border colours.")
    return(NULL)
    } else if(is.null(names(cat.values))||!is.vector(cat.values)){
      message("<< Add sample category >>")
      message("The new category must be a vector of values identifiable with sample names.")
      return(NULL)
    } else {
      # create a dataframe of the new sample categories 
      if (samp.type=="tumour"){
        samp.inf = samp.inf[grep("_T", samp.inf$Sample_Name),]
      } else {
        samp.inf = samp.inf
      }
      rownames(samp.inf) = samp.inf$Sample_Name
      if (length(which(rownames(samp.inf)%in%names(cat.values)))==nrow(samp.inf)){
        if (is.null(cat.levels)){
          cat.levels = c(unique(cat.values[order(cat.values)]))
        }
        new.samp.inf = dplyr::mutate(samp.inf, new_cat = sapply(samp.inf$Sample_Name, function(x) cat.values[x]))
        new.samp.inf[,cat.name] = new.samp.inf$new_cat
        rownames(new.samp.inf) = new.samp.inf$Sample_Name
        new.samp.inf = new.samp.inf[names(cat.values),]
        new.samp.cat=unique(c(samp.cat, cat.name))
        names(new.samp.cat) = c(new.samp.cat)
        samp.lev[[cat.name]] = cat.levels
        cat.cols[[cat.name]] = col[1:length(cat.levels)]
        cat.col.borders[[cat.name]] = col.bord[1:length(cat.levels)]
        names(cat.cols[[cat.name]]) = cat.levels
        names(cat.col.borders[[cat.name]]) = cat.levels
        new.samp.details = list(new.samp.inf, new.samp.cat, samp.lev, cat.cols, cat.col.borders, samp.clust)
        names(new.samp.details) = c("Sample_info", "sample_categories", "sample_levels", "colour_cat", "colour_cat_border", "Sample_Clusters")
        return(new.samp.details)
      } else {
        message("<< Add sample category >>")
        message("The new category vector of values must have the same names as the sample names and be of the same length.")
        return(NULL)
      }
    }
}

# **********
# function to transform a pathways matrix to a pathways list 
# **********
matrix.to.list = function(pws){
  pws.list = list()
  for (pw in colnames(pws)){
    pws.list[[pw]] = rownames(pws)[as.logical(pws[,pw])]
  }
  return(pws.list)
}

# **********
# function to filter out the pathways from a set that are under-represented in your gene list
# **********
prepare.gmt = function(gmt.file, genes.in.data, path.type){
  gmt.pth = gmtPathways(gmt.file)
  hidden = unique(unlist(gmt.pth))
  mat = matrix(NA, dimnames = list(hidden, names(gmt.pth)),
               nrow = length(hidden), ncol = length(gmt.pth))
  for (i in 1:dim(mat)[2]){
    mat[,i] = as.numeric(hidden %in% gmt.pth[[i]])
  }
  hidden1 = intersect(genes.in.data, hidden)
  mat = mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>ifelse(path.type=="Hallmark", 10, 10))]]
  final.list = matrix.to.list(mat)
  return(final.list)
}

# **********
# function to perform a pathway analysis on the genes 
# **********
perform.pathway.analysis = function(genes = NULL,
                                    delta.values = NULL,
                                    fc.values = NULL,
                                    pathway.type = NULL,
                                    subgroup.comp = c("Group1", "Group2"),
                                    gene.type = c("promoter", "All"),
                                    limma.analysis = TRUE,
                                    use.mean = TRUE,
                                    gene.set.filter = FALSE,
                                    min.set=10,
                                    results.dir
                                    ){
  if (pathway.type%in%c("GO", "Reactome")){
    if (!is.null(genes)) {
      entrezids = AnnotationDbi::select(org.Hs.eg.db, keys=genes, columns='ENTREZID', keytype='SYMBOL')
      if (length(which(duplicated(entrezids$SYMBOL)==TRUE))>0){
        entrezids = entrezids[-which(duplicated(entrezids$SYMBOL)==TRUE),]
      }
      if (length(which(is.na(entrezids$ENTREZID)==TRUE))>0){
        entrezids = entrezids[-which(is.na(entrezids$ENTREZID)),]
      }
    }
    if (!is.null(delta.values)){
      entrezids %<>% dplyr::mutate(delta.values = sapply(SYMBOL, function(x) delta.values[x]))
      delta.values = delta.values[which(names(delta.values)%in%entrezids$SYMBOL)]
    }
    if (!is.null(fc.values)){
      entrezids %<>% dplyr::mutate(fc.values = sapply(SYMBOL, function(x) fc.values[x]))
      gene.list = fc.values[which(names(fc.values)%in%entrezids$SYMBOL)]
    }
  }
  enrich.list = list()
  gse.list = list()
  subgroup.title = paste(subgroup.comp, collapse="_vs_")
  if (pathway.type=="GO"){
    enrich.types = c("CC", "MF", "BP")
    for (e.type in enrich.types){
      enrich.list[[paste(pathway.type, e.type, sep="_")]] = enrichGO(gene=entrezids$ENTREZID, OrgDb = org.Hs.eg.db, ont = e.type, pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
      # save the top pathway details to file
      file.name = paste(study_prefix, subgroup.title, ifelse(limma.analysis, "limma", ""), ifelse(use.mean, "mean", "median"), "no_filter",
                        paste0(toupper(substr(gene.type,1,1)), substr(gene.type, 2, nchar(gene.type))), "Regions", paste(pathway.type, e.type, sep="_"), "Top_Pathways.csv", sep="_")
      promoter.pathways = as.data.frame(enrich.list[[paste(pathway.type, e.type, sep="_")]]@result)
      if (!dir.exists(results.dir)){
        create.dir(results.dir)
      }
      write_csv(promoter.pathways, file=paste0(results.dir, file.name))   
    }
    gene.list = gene.list[order(gene.list, decreasing=TRUE)]
    gse.list[[pathway.type]] = gseGO(geneList = gene.list,
                                     ont = "ALL",
                                     keyType = "SYMBOL",
                                     minGSSize = min.set,
                                     maxGSSize = 500,
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "fdr",
                                     OrgDb = org.Hs.eg.db)
    # save the top pathway details to file
    file.name = paste(study_prefix, subgroup.title, ifelse(limma.analysis, "limma", ""), ifelse(use.mean, "mean", "median"), "no_filter",
                      paste0(toupper(substr(gene.type,1,1)), substr(gene.type, 2, nchar(gene.type))), "Regions", pathway.type, "GSEA.csv", sep="_")
    gse.pathways = as.data.frame(gse.list[[pathway.type]]@result)
    if (!dir.exists(results.dir)){
      create.dir(results.dir)
    }
    write_csv(gse.pathways, file=paste0(results.dir, file.name))   
  } else if (pathway.type=="Reactome") {
    enrich.list[[pathway.type]] = enrichPathway(gene=entrezids$ENTREZID, readable=TRUE, maxGSSize = 500, minGSSize = min.set) #nrow(entrezids) 
    # save the top pathway details to file
    file.name = paste(study_prefix, subgroup.title, ifelse(limma.analysis, "limma", ""), ifelse(use.mean, "mean", "median"), "no_filter",
                      paste0(toupper(substr(gene.type,1,1)), substr(gene.type, 2, nchar(gene.type))), "Regions", pathway.type, "Top_Pathways.csv", sep="_")
    promoter.pathways = as.data.frame(enrich.list[[pathway.type]]@result)
    if (!dir.exists(results.dir)){
      create.dir(results.dir)
    }
    write_csv(promoter.pathways, file=paste0(results.dir, file.name)) 
  } else {
    if (pathway.type=="Hallmark"){
      gene.set = paste0(working_dir, "Sample_Sheets/Pathway_Gene_Sets/h.all.v2024.1.Hs.symbols.gmt")
      gene.pth = gmtPathways(gene.set)
      if (gene.set.filter){
        gene.filters = as.data.frame(fread(file=paste0(working_dir, "Sample_Sheets/Pathway_Gene_Sets/Gene_Sets_of_Interest.csv")))
        gene.pth = gene.pth[which(names(gene.pth)%in%gene.filters$Hallmark_Targets)]
      }
    } else if (pathway.type =="C2"){
      if (gene.set.filter){
        gene.set = paste0(working_dir, "Sample_Sheets/Pathway_Gene_Sets/c2.cgp.v2025.1.Hs.symbols.gmt")
        gene.pth = prepare.gmt(gmt.file=gene.set, genes.in.data=names(fc.values), path.type=pathway.type)
        if (gene.set.filter){
          gene.filters = as.data.frame(fread(file=paste0(working_dir, "Sample_Sheets/Pathway_Gene_Sets/Gene_Sets_of_Interest.csv")))
          gene.pth = gene.pth[which(names(gene.pth)%in%gene.filters$CGP_Targets)]
        }
      } else {
        gene.set = paste0(working_dir, "Sample_Sheets/Pathway_Gene_Sets/c2.all.v2024.1.Hs.symbols.gmt")
        gene.pth = prepare.gmt(gmt.file=gene.set, genes.in.data=names(fc.values), path.type=pathway.type)
      }
    } else if (length(grep("Nabet", pathway.type))>0){
      gene.set = paste0(working_dir, "Sample_Sheets/Pathway_Gene_Sets/Nabet_2024_GSEA.gmt")
      gene.pth = gmtPathways(gene.set)
      for(g.pth in names(gene.pth)){
        gene.pth[[g.pth]] = gsub(" ", "", gene.pth[[g.pth]])
      }
      if (length((grep("Hall", pathway.type))>0)){
        gene.pth = gene.pth[grep("Hallmark", names(gene.pth))]
      } else {
        gene.pth = gene.pth[grep("Hallmark", names(gene.pth), invert=TRUE)]
      }
    }
    gene.list = fc.values[order(fc.values, decreasing=TRUE)]
    fgsea.result = fgseaMultilevel(pathways = gene.pth,  #fgsea
                         stats = gene.list,
                         nPermSimple = 10000,
                         scoreType = 'std',
                         minSize = min.set, 
                         maxSize = 500
                         )
    fgsea.result = as.data.frame(fgsea.result[,1:7]) %>% mutate(leadingEdge = unlist(lapply(fgsea.result$leadingEdge, function(x) paste(x, collapse="/"))))
    fgsea.result %<>% dplyr::mutate(ONTOLOGY=rep(ifelse(pathway.type=="Hallmark", ifelse(gene.set.filter, "H_filt", "H"), ifelse(pathway.type=="Nabet_Other", "NO", ifelse(pathway.type=="Nabet_Hall", "NH", ifelse(gene.set.filter, "C2_CGP_filt", "C2")))), nrow(fgsea.result)),
                                    Description = fgsea.result$pathway,
                                    ) %>%
      dplyr::rename(ID=pathway, 
                    setSize=size, 
                    enrichmentScore=ES,
                    pvalue=pval, 
                    p.adjust=padj, 
                    core_enrichment=leadingEdge) %>% 
      dplyr::select(ONTOLOGY, ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, log2err, core_enrichment) %>% dplyr::filter(p.adjust<ifelse(pathway.type=="Hallmark", 0.3, ifelse(length(grep("Nabet", pathway.type))>0, 0.3, 0.1))) %>% dplyr::arrange(p.adjust)
    rownames(fgsea.result) = fgsea.result$ID
    gse.list[[pathway.type]] = new("gseaResult", 
                                   result=fgsea.result,
                                   organism = "Homo sapien",
                                   setType = pathway.type,
                                   geneSets = gene.pth,
                                   geneList = gene.list,
                                   keytype="SYMBOL",
                                   permScores = matrix(,nrow=0, ncol=0),
                                   params = list(pvalueCutoff=ifelse(pathway.type=="Hallmark", 0.3, ifelse(length(grep("Nabet", pathway.type))>0, 0.3, 0.1)), eps=1e-10, pAdjustMethod="BH", exponent=1, minGSSize=min.set, MaxGSSize=500), 
                                   gene2Symbol = "",
                                   readable= TRUE,
                                   termsim = matrix(,nrow=0, ncol=0),
                                   method = "",
                                   dr = list())
    # save the top pathway details to file
    file.name = paste(study_prefix, subgroup.title, ifelse(limma.analysis, "limma", ""), ifelse(use.mean, "mean", "median"), "no_filter",  
                      paste0(toupper(substr(gene.type,1,1)), substr(gene.type, 2, nchar(gene.type))), "Regions", ifelse(gene.set.filter, paste0("filtered_", pathway.type), pathway.type), "GSEA.csv", sep="_")
    gse.pathways = as.data.frame(gse.list[[pathway.type]]@result)
    if (!dir.exists(results.dir)){
      create.dir(results.dir)
    }
    write_csv(gse.pathways, file=paste0(results.dir, file.name))   
  }
}

# ***************
# function to draw just the heatmap annotations 
# ***************
draw.clinical.annotations = function(samp.info,
                                     samp.cat,
                                     samp.lev,
                                     cat.col,
                                     matrix.comparison = NULL,
                                     hm.matrix = NULL,
                                     hm.matrix.col = colorRamp2(c(-0.6,0,0.6), c("blue", "white", "red")),
                                     hm.lab = NULL,
                                     hm.r.pad = 13,
                                     extra.matrix = NULL,
                                     extra.matrix.col = colorRamp2(c(0,1), c("white", "dark blue")),
                                     extra.lab = NULL,
                                     extra.r.pad = 13,
                                     push.list = NULL
                                     ){
  annot.font.size = 8
  samp.font.size = 5
  row.font.size = 7
  title.font.size = 14
  legend.title.cex = 1.1 
  legend.text.cex = 1 
  point.type = 15
  point.cex = 2.3 
  
  hm.zero = matrix(nrow=0, ncol = nrow(samp.info))
  
  colnames(hm.zero) = samp.info$Sample_Name
  annotations = samp.info[which(samp.info$Sample_Name %in% colnames(hm.zero)), samp.cat]
  rownames(annotations) = samp.info$Sample_Name
  anno.cols = list()
  for (cat in samp.cat){
    samp.cols = cat.col[[cat]]
    samp.cols = samp.cols[samp.lev[[cat]]]
    anno.cols[[cat]] = c(samp.cols)
  }
  names(anno.cols) = colnames(annotations)
  anno.cols$Meth_Cluster = anno.cols$Meth_Cluster[c("Group1", "Group2", "Group3", "Group4")]
  par(mfrow=c(3,1))
  plot.new()
  
  # first plot annotations
  left.pad = 30
  right.pad = 20.4
  par(mai=c(0.1, 0.1, 0.1, 0.1))
  ha = HeatmapAnnotation(df = annotations,
                         col=anno.cols, 
                         show_legend=FALSE,
                         show_annotation_name = TRUE,
                         annotation_name_side = "left",
                         annotation_name_gp = gpar(fontsize=annot.font.size)#,
                         )

  ht.annot = Heatmap(hm.zero, 
                     cluster_rows=FALSE,
                     cluster_columns=FALSE,
                     row_names_side="left",
                     row_names_gp=gpar(fontsize=row.font.size),
                     bottom_annotation=ha, 
                     show_column_names = FALSE, 
                     column_names_gp = gpar(fontsize=samp.font.size)
                     )
  annot.grob = grid.grabExpr(draw(ht.annot, padding = unit(c(1,left.pad,0,right.pad),"mm"))) 
  
  if (!is.null(hm.matrix)){
    ht = Heatmap(hm.matrix, 
                 col = hm.matrix.col,
                 cluster_rows=TRUE,
                 cluster_columns=FALSE,
                 show_row_names=ifelse(nrow(hm.matrix)<50, TRUE, FALSE),
                 row_names_side="left",
                 row_names_gp=gpar(fontsize=row.font.size),
                 show_column_names = FALSE,
                 heatmap_legend_param = list(title=hm.lab, 
                                             title_gp=gpar(fontsize=6), 
                                             labels_gp=gpar(fontsize=5), 
                                             legend_height=unit(2,"cm"), 
                                             grid_width=unit(2,"mm"),
                                             just=c("right", "top"))
                 )
    if (nrow(hm.matrix)<50) {
      l.pad = left.pad - (max(nchar(rownames(hm.matrix)))+1.9)
    } else {
      l.pad = left.pad
    }
    ht.grob = grid.grabExpr(draw(ht, padding = unit(c(1,l.pad,2,hm.r.pad),"mm")))
  } else {
    ht.grob = NULL
  }
    
  if (!is.null(extra.matrix)){
    ht.extra = Heatmap(extra.matrix, 
                       col = extra.matrix.col,
                       cluster_rows=TRUE,
                       cluster_columns=FALSE,
                       show_row_names=ifelse(nrow(extra.matrix)<50, TRUE, FALSE),
                       row_names_side="left",
                       row_names_gp=gpar(fontsize=row.font.size),
                       column_title=paste0(study_name, " - Clinical Data ", matrix.comparison),
                       column_title_gp = gpar(fontsize=title.font.size),
                       show_column_names = FALSE,
                       heatmap_legend_param = list(title=extra.lab, 
                                                   title_gp=gpar(fontsize=6), 
                                                   labels_gp=gpar(fontsize=5), 
                                                   legend_height=unit(2,"cm"), 
                                                   grid_width=unit(2,"mm"),
                                                   just=c("right", "top"))
                       )
    if (nrow(extra.matrix)<50) {
      l.pad = left.pad - (max(nchar(rownames(extra.matrix)))+1.9)
    }
    
    extra.grob = grid.grabExpr(draw(ht.extra, padding = unit(c(1,l.pad,1,extra.r.pad),"mm")))
  } else {
    extra.grob = NULL
  }
  return(list(annot=annot.grob, ht=ht.grob, extra=extra.grob))
}

# ***************
# Create function to extract the probe information related to the gene of interest
# ***************
load.gene.info = function(probe.info, 
                          gene,
                          levels.of.interest = c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR")){
  gene.data = NULL
  columns = c("Name", "Genome_Build", "CHR", "MAPINFO", "Strand", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "UCSC_RefGene_Accession")
  gene.data = probe.info[c(grep(gene, probe.info$UCSC_RefGene_Name, fixed=TRUE)), c(columns)]
  colnames(gene.data) = c("TargetID", "RefGene", "Chromosome", "Position", "Strand", "Hugo_Symbol", "Annotation", "Accession_Number")
  gene.data$Annotation = as.character(gene.data$Annotation)
  gene.rows = NULL
  for (i in 1:nrow(gene.data)){
    if (gene %in% unlist(strsplit(as.character(gene.data$Hugo_Symbol[i]),";", fixed=TRUE))){
      gene.rows = c(gene.rows, i)
    }
  }
  gene.data = gene.data[gene.rows,]
  gene.data$Hugo_Symbol = as.character(gene.data$Hugo_Symbol)
  accessions = table(unlist(strsplit(as.character(gene.data$Accession_Number), ";", fixed=TRUE)))
  accession.number = names(accessions)[order(accessions, decreasing=TRUE)]
  rec = 0
  for (k in 1:nrow(gene.data)){
    rec_accessions = unlist(strsplit(as.character(gene.data$Accession_Number[k]),";", fixed=TRUE))
    if (length(rec_accessions)==1){
      access.num = 1
    } else {
      if (length(which(rec_accessions==accession.number[1]))>0){
        if(unique(unlist(strsplit(as.character(gene.data$Hugo_Symbol[k]),";", 
                                  fixed=TRUE))[which(rec_accessions==accession.number[1])])==gene){
          access.num = which(rec_accessions==accession.number[1])
        } else if (unique(unlist(strsplit(as.character(gene.data$Hugo_Symbol[k]),";", 
                                          fixed=TRUE))[which(rec_accessions==accession.number[2])])==gene){
          access.num = which(rec_accessions==accession.number[2])
        } else {
          access.num = which(rec_accessions==accession.number[3])
        }
      } else if(length(which(rec_accessions==accession.number[2]))>0) {
        if(unique(unlist(strsplit(as.character(gene.data$Hugo_Symbol[k]),";", 
                                  fixed=TRUE))[which(rec_accessions==accession.number[2])])==gene){
          access.num = which(rec_accessions==accession.number[2])
        } else {
          access.num = which(rec_accessions==accession.number[3])
        }
      }
    }
    gene.num = which(unlist(strsplit(as.character(gene.data$Hugo_Symbol[k]),";", fixed=TRUE))==gene)
    gene.data$Hugo_Symbol[k] = unique(unlist(strsplit(as.character(gene.data$Hugo_Symbol[k]),";", fixed=TRUE))[intersect(access.num, gene.num)])
    
    if (length(unique(unlist(strsplit(as.character(gene.data$Annotation[k]), ";", fixed = TRUE))[intersect(access.num, gene.num)]))>1){
      if (k==1){
        gene.data$Annotation[k] = unique(unlist(strsplit(as.character(gene.data$Annotation[k]), ";", fixed = TRUE))[intersect(access.num, gene.num)])[1]
      } else {
        gene.data$Annotation[k] = gene.data$Annotation[k-1]
      }
    } else {
      gene.data$Annotation[k] = unique(unlist(strsplit(as.character(gene.data$Annotation[k]), ";", fixed = TRUE))[intersect(access.num, gene.num)])
    }
    gene.data$Annotation[k] = ifelse(gene.data$Annotation[k]=="ExonBnd", "Body", gene.data$Annotation[k])
  }
  gene.data$PosLabel = vector(mode="character", length=nrow(gene.data))
  
  for (i in 1:nrow(gene.data)){
    gene.data$PosLabel[i] = paste0("Chr", gene.data$Chromosome[i], ":", gene.data$Position[i])
  }
  
  for (j in 1:ncol(gene.data)){
    if (is.factor(gene.data[,j])) {
      gene.data[,j] = droplevels(gene.data[,j])
    }
  }
  gene.data$Position = as.integer(gene.data$Position)
  if (nrow(gene.data)>0){
    if (as.character(gene.data$Annotation[order(gene.data$Position, decreasing=T)])[1] %in% c("TSS1500", "TSS200", "5'UTR")){
      order.gene = TRUE
    } else{
      order.gene = FALSE
    }
  }
  gene.data$PosLabel = factor(gene.data$PosLabel, levels=gene.data$PosLabel[order(gene.data$Position, decreasing=order.gene)])
  gene.data = gene.data[order(match(gene.data$PosLabel, levels(gene.data$PosLabel))),]
  gene.data %<>% dplyr::filter(Annotation %in% levels.of.interest)
  gene.data$Annotation = factor(gene.data$Annotation, levels=c(levels.of.interest))
  gene.data = gene.data[order(gene.data$Annotation),]
  return(gene.data)
}


# ***************
# Function to create heatmap of beta values for samples with RNASeq data for the gene of interest
# ***************
draw.gene.heatmap = function(meth.data = meth_norm,
                             samp.info = sample_info,
                             RNA.info = NULL, 
                             zscore=FALSE,
                             LoH.info = NULL, 
                             cat.cols = colour_cat,
                             gene.info = gene_info,
                             data.type = NULL,
                             gene = NULL,
                             meth.sub = "",
                             gene.meth.sub = "Meth_Cluster",
                             status = NULL,
                             add.annos = "Tumour_Content"){
  if(!is.null(gene)){
    if (is.null(gene.info) && !gene=="all"){
      gene.info=list()
      gene.info[[gene]] = load.gene.info(Probe_info, gene)
    }
    if (gene=="all" && !is.null(gene.info)) {
      gene.info[["all"]] = NULL
      for (gene.name in names(gene.info)){
        gene.info[["all"]] = rbind(gene.info[["all"]], gene.info[[gene.name]])
      }
    }
    gene.samp = samp.info
    gene.locus.name = colnames(gene.samp)[grep(paste0("Gene_Locus_", gene), colnames(gene.samp))]
    if (length(gene.locus.name)>0){
      gene.samp %<>% dplyr::select(Sample_Name, Sample_Group, gene.locus.name) # , Sample_Cluster
    } else if (gene=="all") {
      if (!gene.meth.sub==""){
        gene.samp %<>% dplyr::mutate(Gene_Meth_Subgroup = samp.info[,gene.meth.sub])
        if (!meth.sub[1]==""){
          num.sub = 1
          for (m.sub in meth.sub){
            gene.samp %<>% dplyr::mutate(Methylation_Subgroup=samp.info[,meth.sub[num.sub]])
            colnames(gene.samp)[which(colnames(gene.samp)=="Methylation_Subgroup")] = paste("Methylation_Subgroup", num.sub)
            num.sub = num.sub+1
          }
          gene.samp %<>% dplyr::select(Sample_Name, Sample_Group, Gene_Meth_Subgroup, all_of(colnames(gene.samp)[grep("Methylation_Subgroup", colnames(gene.samp))]), all_of(add.annos))
        } else {
          gene.samp %<>% dplyr::select(Sample_Name, Sample_Group, Gene_Meth_Subgroup, all_of(add.annos))
        }
      } else {
        if (!meth.sub[1]==""){
          num.sub = 1
          for (m.sub in meth.sub){
            gene.samp %<>% dplyr::mutate(Methylation_Subgroup=samp.info[,meth.sub[num.sub]])
            colnames(gene.samp)[which(colnames(gene.samp)=="Methylation_Subgroup")] = paste("Methylation_Subgroup", num.sub)
            num.sub = num.sub+1
          }
          gene.samp %<>% dplyr::select(Sample_Name, Sample_Group, all_of(colnames(gene.samp)[grep("Methylation_Subgroup", colnames(gene.samp))]), all_of(add.annos))
        } else {
          gene.samp %<>% dplyr::select(Sample_Name, Sample_Group, all_of(add.annos))
        }
      }
    } else {
      if (!meth.sub[1]==""){
        num.sub = 1
        for (m.sub in meth.sub){
          gene.samp %<>% dplyr::mutate(Methylation_Subgroup=samp.info[,meth.sub[num.sub]])
          colnames(gene.samp)[which(colnames(gene.samp)=="Methylation_Subgroup")] = paste("Methylation_Subgroup", num.sub)
          num.sub = num.sub+1
        }
        gene.samp %<>% dplyr::select(Sample_Name, Sample_Group, all_of(colnames(gene.samp)[grep("Methylation_Subgroup", colnames(gene.samp))]), all_of(add.annos))
      } else {
        gene.samp %<>% dplyr::select(Sample_Name, Sample_Group, all_of(add.annos))
      }
    }
    gene.samp$Sample_Name = as.character(gene.samp$Sample_Name)
    gene.meth = meth.data[which(rownames(meth.data)%in%c(gene.info[[gene]]$TargetID)), gene.samp$Sample_Name]
    probe.order = as.character(gene.info[[gene]]$TargetID[which(gene.info[[gene]]$TargetID%in%rownames(gene.meth))])
    
    if (!is.matrix(gene.meth)){
      gene.meth = as.matrix(t(gene.meth), nrow=1, byrow=TRUE)
      rownames(gene.meth) = probe.order
    } else {
      gene.meth = gene.meth[probe.order,]
    }
    if (!is.null(LoH.info)){
      LoH.info = LoH.info[,gene.samp$Sample_Name]
      Chr.cols = ifelse(unique(gene.info[[gene]]$Chromosome) %in% c("X","Y"), "",
                        ifelse(gene.info[[gene]]$Position[1] < 
                                 centromere_info$Centromere_start[which(centromere_info$CHR==gene.info[[gene]]$Chromosome[1])],
                               paste0("LoH_Chr", unique(gene.info[[gene]]$Chromosome)),
                               paste0("LoH_Chr", unique(gene.info[[gene]]$Chromosome), ".5")))
      rownames(LoH.info) = paste0("LoH_Chr", rownames(LoH.info))
      for (col in Chr.cols){
        if (col %in% rownames(LoH.info)){
          gene.samp[,col] = vector(mode="integer", length=nrow(gene.samp))
          gene.samp[,col] = unlist(LoH.info[which(rownames(LoH.info)==col),])
        }
      }
    } else {
      Chr.cols = "Not_Used"
    }
    if(!is.null(RNA.info)){
      if (!gene=="all"){
        RNA.gene = gene
      } else {
        if (!is.null(status)){
          if (status=="parpi"){
            RNA.gene = c(names(gene.info), "HUWE1")
          } else {
            RNA.gene = names(gene.info)
          }
        } else {
          RNA.gene = names(gene.info)
        }
        RNA.gene = RNA.gene[grep("all", RNA.gene, invert = TRUE)]
      }
      if (length(RNA.gene)==1){
        RNA.info %<>% dplyr::filter(`Associated Gene Name`== RNA.gene)
      } else {
        RNA.info %<>% dplyr::filter(`Associated Gene Name`%in% RNA.gene)
      }
      for (gene.name in RNA.gene){
        RNA.samp.info = data.frame(Sample_Name = colnames(RNA.info[grep("_T", colnames(RNA.info))]),
                                   RNA_Exp=as.vector(as.matrix(RNA.info[which(RNA.info$`Associated Gene Name`==gene.name), grep("_T", colnames(RNA.info))], ncol=1)))
        RNA.samp.info %<>% dplyr::filter(Sample_Name%in%colnames(gene.meth))
        if (zscore){
          RNA.samp.info$RNA_zScore = (RNA.samp.info$RNA_Exp-(mean(RNA.samp.info$RNA_Exp)))/sd(RNA.samp.info$RNA_Exp)
        } else {
          RNA.samp.info$RNA_Log2 = log2(RNA.samp.info$RNA_Exp+1)
        }
        gene.samp = merge(gene.samp, RNA.samp.info, by.x="Sample_Name", by.y="Sample_Name", all.x=TRUE)
        if (zscore){
          gene.samp$RNA_zScore[which(is.infinite(gene.samp$RNA_zScore)==TRUE)] = 9*-1
          colnames(gene.samp)[which(colnames(gene.samp)=="RNA_zScore")] = paste(gene.name, "RNA_zScore", sep="_")
        } else {
          colnames(gene.samp)[which(colnames(gene.samp)=="RNA_Log2")] = paste(gene.name, "RNA_Log2", sep="_")
        }
        colnames(gene.samp)[which(colnames(gene.samp)=="RNA_Exp")] = paste(gene.name, "RNA_Exp", sep="_")
      }
    }
    rownames(gene.samp) = gene.samp$Sample_Name
    gene.samp = gene.samp[colnames(gene.meth),]
    if (!is.null(RNA.info)){
      RNA.values = which(!is.na(gene.samp[,c(paste(RNA.gene[1], ifelse(zscore, "RNA_zScore", "RNA_Log2"), sep="_"))])) # 
    }
    if (!is.null(RNA.info) && !gene=="all"){
      correlations = apply(gene.meth[,gene.samp$Sample_Name[RNA.values]], 
                           1, function(x) {cor.test(x, gene.samp[RNA.values, paste(gene.name, "RNA_Exp", sep="_")], 
                                                    alternative="two.sided", method="pearson")})
    } else {
      correlations = NULL
    }
    if (status %in% c("remove")){ # c("oncogene", "treatment")
      # read in the copy number for the oncogenes from 77 WGS SCLC analysis
      load(WGS_VEP_data)
      CNV_data %<>% dplyr::filter(Hugo_Symbol%in%gsub("MYCL1", "MYCL", names(gene.info)))
      study.data = as.data.frame(fread(file=analysis_file))
      study.data %<>% dplyr::filter(donorLabel%in%gsub("_T", "", samp.info$Sample_Name[which(samp.info$RNA=="YES")]), !ascatAnalysis=="")
      replicas = unique(study.data$donorLabel[which(duplicated(study.data$donorLabel))])
      keep = c()
      for (donor in replicas){
        keep = c(keep, which(study.data$donorLabel==donor&study.data$testSequencingPlatform=="NovaSeq"))
      }
      study.data = rbind(dplyr::filter(study.data, !donorLabel%in%replicas), study.data[keep,])
      study.data %<>% dplyr::mutate(donorLabel=paste0(donorLabel, "_T"))
      for (gene.name in gsub("MYCL1", "MYCL", names(gene.info))){
        gene.samp %<>% dplyr::mutate(CN_gene = sapply(Sample_Name, function(x) ifelse(nrow(dplyr::filter(CNV_data, paste(DonorLabel, "T", sep = "_")==x, Hugo_Symbol==gene.name))>0,
                                                                                      ifelse(dplyr::filter(CNV_data, paste(DonorLabel, "T", sep = "_")==x, Hugo_Symbol==gene.name)$Copy_Number_Variant%in%c("Del", "cnLOH"),
                                                                                             "Loss/LOH",
                                                                                             ifelse(dplyr::filter(CNV_data, paste(DonorLabel, "T", sep = "_")==x, Hugo_Symbol==gene.name)$Copy_Number_Variant%in%c("Gain", "HighGain"),
                                                                                                    ifelse(dplyr::filter(study.data, donorLabel==x)$ascatPloidyScore>=2.7,
                                                                                                           ifelse(dplyr::filter(CNV_data, paste(DonorLabel, "T", sep = "_")==x, Hugo_Symbol==gene.name)$Copynumber>=9, 
                                                                                                                  "HighGain",
                                                                                                                  ifelse(dplyr::filter(CNV_data, paste(DonorLabel, "T", sep = "_")==x, Hugo_Symbol==gene.name)$Copynumber%in%c(5,6,7,8),
                                                                                                                         "Gain",
                                                                                                                         "Het")),
                                                                                                           dplyr::filter(CNV_data, paste(DonorLabel, "T", sep = "_")==x, Hugo_Symbol==gene.name)$Copy_Number_Variant),
                                                                                                    "Het")), 
                                                                                      "CN_NA")))
        colnames(gene.samp)[which(colnames(gene.samp)=="CN_gene")] = paste("CN", gene.name, sep="_")
      }
    }
    if (!is.null(LoH.info)&&!is.null(RNA.info)&&length(RNA.gene)==1){
      annotation.order = c("RNA_Exp", Chr.cols, gene.locus.name)
    } else if (!is.null(LoH.info)){
      annotation.order = c(Chr.cols, gene.locus.name) 
    } else if (!is.null(RNA.info)&&length(RNA.gene)==0){
      annotation.order = c(colnames(gene.samp)[grep(ifelse(zscore, "RNA_zScore", "RNA_Log2"), colnames(gene.samp))], add.annos, colnames(gene.samp)[grep("Methylation_Subgroup", colnames(gene.samp))]) 
    } else if (!is.null(RNA.info)&&length(RNA.gene)>=1){
      if (status %in% c("remove")){ # c("oncogene", "treatment")
        if (gene.meth.sub==""){
          annotation.order = c(colnames(gene.samp)[grep(ifelse(zscore, "RNA_zScore", "RNA_Log2"), colnames(gene.samp))], colnames(gene.samp)[grep("CN_MYC|CN_DLL|CN_PLCG|CN_BCL|CN_SLFN|CN_KDM|CN_SEZ|CN_CD276|CN_TACSTD", colnames(gene.samp))], add.annos, colnames(gene.samp)[grep("Methylation_Subgroup", colnames(gene.samp))]) 
        } else {
          annotation.order = c(colnames(gene.samp)[grep(ifelse(zscore, "RNA_zScore", "RNA_Log2"), colnames(gene.samp))], colnames(gene.samp)[grep("CN_MYC|CN_DLL|CN_PLCG|CN_BCL", colnames(gene.samp))], add.annos, colnames(gene.samp)[grep("Methylation_Subgroup", colnames(gene.samp))], "Gene_Meth_Subgroup")
        }
      } else {
        if (gene.meth.sub==""){
          annotation.order = c(colnames(gene.samp)[grep(ifelse(zscore, "RNA_zScore", "RNA_Log2"), colnames(gene.samp))], add.annos, colnames(gene.samp)[grep("Methylation_Subgroup", colnames(gene.samp))]) 
        } else {
          annotation.order = c(colnames(gene.samp)[grep(ifelse(zscore, "RNA_zScore", "RNA_Log2"), colnames(gene.samp))], add.annos, colnames(gene.samp)[grep("Methylation_Subgroup", colnames(gene.samp))], "Gene_Meth_Subgroup") 
        }
      }
    } else {
      annotation.order = c("Sample_Group", add.annos, colnames(gene.samp)[grep("Methylation_Subgroup", colnames(gene.samp))]) 
    }
    annotations = dplyr::select(gene.samp, which(colnames(gene.samp) %in% annotation.order))
    annotation.order = annotation.order[which(annotation.order%in%colnames(annotations))]
    annotations = annotations[,annotation.order]
    for (col.name in colnames(annotations)[grep(ifelse(zscore, "RNA_zScore", "RNA_Log2"), colnames(annotations))]){ 
      annotations[which(is.na(annotations[,col.name])), col.name] = 0
    }
    colnames(annotations)[grep("Methylation_Subgroup", colnames(annotations))] = meth.sub
    if (!gene.meth.sub==""){
      colnames(annotations)[which(colnames(annotations)=="Gene_Meth_Subgroup")] = gene.meth.sub
    }
    anno.cols = vector("list", ncol(annotations))
    names(anno.cols) = colnames(annotations)
    if (zscore){
      rna.cols = c("dodgerblue4", "white", "firebrick")  
    } else {
      rna.cols = c("navy", "darkcyan", "gold")
    }
    legend.labels = list()
    for (variable in names(anno.cols)){
      if (length(grep(ifelse(zscore, "RNA_zScore", "RNA_Log2"),variable))>0){ 
        anno.cols[[variable]] = colorRamp2(c(ifelse(zscore, ifelse(min(annotations[,variable])<=0, min(annotations[,variable]), (2*-1)), 0), 
                                             ifelse(zscore, 0, 4), 
                                             ifelse(zscore, ceiling(max(annotations[RNA.values,variable])), 8)),
                                           rna.cols)
      } else if (variable==Chr.cols){
        anno.cols[[Chr.cols]] = colorRamp2(c(0,1), c("white", "dark blue"))
      } else if (length(grep("CN_",variable))>0){ #MYC|CN_DLL|CN_PLCG2|CN_BCL
        anno.cols[[variable]] = c(Del="#0c2860", cnLOH="#8ed1db", `Loss/LOH`="#8ed1db", Het="#f0eadc", Gain="#fba194", HighGain="#bc0f0b", CN_NA="white")
        annotations[,variable] = droplevels(factor(annotations[,variable], levels=names(anno.cols[[variable]])))
      } else {
        anno.cols[[variable]] = cat.cols[[variable]]
        annotations[,variable] = factor(annotations[,variable], levels=names(anno.cols[[variable]]))
      }
    }
    if (!is.null(LoH.info)){
      colnames(annotations)[which(colnames(annotations)==gene.locus.name)] = "Gene_Locus"
      colnames(annotations)[which(colnames(annotations)==Chr.cols)] = paste("Avg", Chr.cols, sep=" ")
    }
    if (!is.null(RNA.info)){
      for (r.gene.name in RNA.gene){
        annotations[,paste(r.gene.name, ifelse(zscore, "RNA_zScore", "RNA_Log2"), sep= "_")] = as.numeric(annotations[,paste(r.gene.name, ifelse(zscore, "RNA_zScore", "RNA_Log2"), sep= "_")]) 
      }
    }
    annot.font.size = 12
    samp.font.size = 8
    title.font.size = 10
    legend.title.cex = 12
    legend.text.cex = 10
    point.type = 15
    point.cex = 1.6 #2.5
    legends.to.show = sapply(colnames(annotations), function(x) ifelse(length(grep("zScore|Exp|Log2|CN_", x))>0, FALSE, TRUE))
    rna.legends = intersect(which(legends.to.show==FALSE), grep("zScore|Exp|Log2", colnames(annotations)))
    legends.to.show[rna.legends[1]]=TRUE
    colnames(annotations) = gsub("_RNA_", " ", colnames(annotations))
    cn.legends = intersect(which(legends.to.show==FALSE), grep("CN_", colnames(annotations)))
    legends.to.show[cn.legends[1]]=TRUE
    colnames(annotations) = gsub("_", " ", colnames(annotations))
    names(anno.cols) = colnames(annotations)
    # change paper region name
    colnames(annotations)[which(colnames(annotations)=="Paper Regions M")] = "Meth Regions Cluster"
    names(anno.cols)[which(names(anno.cols)=="Paper Regions M")] = "Meth Regions Cluster" # "Log2(exp+1)"
    names(legends.to.show)[which(names(legends.to.show)=="Paper Regions M")] = "Meth Regions Cluster" # "Log2(exp+1)"
    legend.param = list(labels_gp = gpar(fontsize=legend.text.cex),
                        title_gp = gpar(fontsize=legend.title.cex, fontface="bold"),
                        grid_height=unit(5,"mm"), 
                        grid_width=unit(5, "mm")) #, at = legend.labels)
    ha = HeatmapAnnotation(df = annotations,
                           col=anno.cols, 
                           show_legend=legends.to.show,
                           show_annotation_name = TRUE,
                           annotation_name_side = "left",
                           annotation_height=unit(ncol(annotations)*4, "mm"),
                           which = "column",
                           annotation_name_gp = gpar(fontsize=annot.font.size),
                           annotation_legend_param = legend.param,
                           na_col = na_colour
                           )
    colourscale = c("#bc0f0b","#f54642", "#f7716e", "#fbbbb1", "#e5d6b3", "#fadeb7")
    names(colourscale)=c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR")
    gene.colourscale = list(Annotation = colourscale)
    gene.info[[gene]] %<>% dplyr::filter(TargetID%in%probe.order)
    gene.annos = data.frame(TargetID = gene.info[[gene]]$TargetID,
                            Annotation = gene.info[[gene]]$Annotation)
    rownames(gene.annos) = gene.annos$TargetID
    gene.annos %<>% dplyr::select(Annotation)
    ha.gene = rowAnnotation(df = gene.annos,
                            col=gene.colourscale, 
                            show_legend=TRUE,
                            show_annotation_name = TRUE,
                            annotation_width = unit(4, "mm"),
                            annotation_name_side = "top",
                            annotation_name_gp = gpar(fontsize=annot.font.size),
                            annotation_legend_param = list(labels_gp = gpar(fontsize=legend.text.cex),
                                                           title_gp = gpar(fontsize=legend.title.cex, fontface="bold"),
                                                           grid_height=unit(5,"mm"), 
                                                           grid_width=unit(5,"mm")),
                            na_col = na_colour
                            )
    row.font.size = ifelse(nrow(gene.meth)>90, 4, 
                           ifelse(nrow(gene.meth)>50, 5, 
                                  ifelse(nrow(gene.meth)>25, 6, annot.font.size)))
    row.font.size = as.integer(row.font.size)
    if(length(RNA.gene)>1){
      row.cut = c()
      num = 0
      for (gene.name in names(gene.info)[grep("all", names(gene.info), invert=TRUE)]){
        num = length(which(gene.info[[gene.name]]$TargetID %in% rownames(gene.meth)))
        row.cut = c(row.cut, rep(gene.name, num))
      }
      row.cut = data.frame(index=c(seq(1,1,length(row.cut))),
                           gene = row.cut)
      row.cut$gene = factor(row.cut$gene, levels = names(gene.info)[grep("all", names(gene.info), invert=TRUE)])
      
    } else {
      row.cut = rep(gene, nrow(gene.meth))
    }
    gene.ht = Heatmap(gene.meth,
                      col=colorRamp2(c(0,0.5,1), c("#082196", "#fefdd2", "#c00e02")),
                      cluster_rows=FALSE,
                      cluster_columns=FALSE,
                      row_names_side="left",
                      row_names_gp=gpar(fontsize=row.font.size),
                      top_annotation=ha, 
                      column_split = gene.samp$Methylation_Subgroup,
                      column_title=paste(study_name, "-", gene, data.type, "Probe Profile."),
                      column_title_gp = gpar(fontsize=title.font.size),
                      show_column_names = FALSE,
                      row_split = row.cut,
                      heatmap_legend_param = list(title="Beta Values", 
                                                  title_gp=gpar(fontsize=legend.title.cex, fontface="bold"), 
                                                  labels_gp=gpar(fontsize=legend.text.cex), 
                                                  legend_height=unit(2,"cm"), 
                                                  grid_width=unit(5,"mm"),
                                                  just=c("right", "top"))
    ) + ha.gene
    gene.grob = grid.grabExpr(draw(gene.ht, gap = unit(1,"mm"), padding = unit(c(1,15,1,10),"mm")))
    gene.ht.info=list(heatmap=gene.grob, meth_info=gene.meth, corrs = correlations)
    return(heatmap=gene.ht.info)
  } else {
    message("<< Plot Gene Matrix. >>")
    message("This function must be supplied with a Gene Name in the format of Hugo symbol or 'all' for all genes in the list to be displayed")
  }
}

# ***************
# Function to create heatmap of beta values for samples with RNASeq data for the gene of interest
# ***************
draw.meth.rna.box.corr.plot = function(meth.data = NULL,
                                       samp.info = sample_info,
                                       RNA.info = NULL, #RNASeq_TMM,
                                       cnv.data = NULL,
                                       col.group = "Meth_Cluster",
                                       cat.cols = colour_cat,
                                       cat.cols.bold = colour_cat_border,
                                       gene.info = NULL,
                                       genes = NULL,
                                       gene.type = "",
                                       violin = FALSE,
                                       plot.sd=0.2,
                                       multi.test = "anova",
                                       single.test = "t.test",
                                       all.comp = FALSE
                                       ){
  if (is.null(meth.data) & is.null(RNA.info)) {
    stop("Must include either Methylation Beta values or RNASeq TMM data.")
  }
  if (!col.group %in% colnames(samp.info)) {
   stop(paste0("Colour group for samples of interest must be included in samp.info as '", col.group, "'"))
  }
  if (!is.null(gene.info)) {
    if (!is.null(meth.data)){
      if (TRUE %in% c(genes %in% names(gene.info))){
        promoter.regions = c("TSS1500", "TSS200", "5`UTR", "1stExon")
        meth.promoter = c()
        meth.body = c()
        meth.sd.promoter = c()
        meth.sd.body = c()
        for (gene in genes) {
          if (gene %in% names(gene.info)){
            gene.probes = gene.info[[gene]]$TargetID
            gene.data = meth.data[gene.probes,]
            gene.promoter.probes = gene.info[[gene]]$TargetID[which(gene.info[[gene]]$Annotation%in%promoter.regions)]
            meth.promoter = c(meth.promoter, gene.promoter.probes)
            gene.body.probes = gene.info[[gene]]$TargetID[which(!gene.info[[gene]]$Annotation%in%promoter.regions)]
            meth.body = c(meth.body, gene.body.probes)
            gene.data = meth.data[gene.info[[gene]]$TargetID,]
            if (length(gene.promoter.probes)>0){
              if (length(gene.promoter.probes)>1){
                meth.sd.promoter = c(meth.sd.promoter, apply(gene.data[gene.promoter.probes,], 1, sd))
              } else {
                meth.sd.promoter = c(meth.sd.promoter, sd(gene.data[gene.promoter.probes,]))
                names(meth.sd.promoter)[length(meth.sd.promoter)] = gene.promoter.probes
              }
            }
            if (length(gene.body.probes)>0){
              if (length(gene.body.probes)>1){
                meth.sd.body = c(meth.sd.body, apply(gene.data[gene.body.probes,], 1, sd))
              } else {
                meth.sd.body = c(meth.sd.body, sd(gene.data[gene.body.probes,]))
                names(meth.sd.body)[length(meth.sd.body)] = gene.body.probes
              }
            }
          }
        }
        no.sd.promoter = ifelse(length(meth.sd.promoter)>0,
                                ifelse(length(which(meth.sd.promoter>=plot.sd))>0, FALSE, TRUE),
                                TRUE)
        no.sd.body = ifelse(length(meth.sd.body)>0, 
                            ifelse(length(which(meth.sd.body>=plot.sd))>0, FALSE, TRUE),
                            TRUE)
        meth.sd = c(meth.sd.promoter, meth.sd.body)
        meth.sd = meth.sd[which(meth.sd>=plot.sd)]
      } else {
        meth.data = NULL
      }
    } else {
      stop("If including gene.info then must supply Methylation Beta values in meth.data.")
    }
  } else {
    if (!is.null(meth.data)){
      stop("In order to preceed with methylation boxplots. Gene/s of interest must be included in gene.info")
    }
  }
  if (!is.null(RNA.info)){
    if (TRUE %in% c(genes %in% RNA.info$`Associated Gene Name`)){
      RNA.info %<>% dplyr::filter(`Associated Gene Name`%in%c(genes))
      rownames(RNA.info) = RNA.info$`Associated Gene Name`
      RNA.info = RNA.info[,colnames(RNA.info)[grep("_T", colnames(RNA.info))]]
      if (!nrow(RNA.info)>0){
        RNA.info=NULL
      }
    } else {
      RNA.info = NULL
    } 
  }
  if (!is.null(meth.data) & !is.null(RNA.info)) {
    gene.avg = data.frame(Gene = rep(gene.type, ncol(meth.data)),
                          Sample_Name = colnames(meth.data),
                          Gp_by = sapply(colnames(meth.data),  function(x) samp.info[,col.group][which(samp.info$Sample_Name==x)]),
                          Meth_Promoter = sapply(colnames(meth.data), function(x) ifelse(length(meth.promoter)>0, 
                                                                                         ifelse(gene.type%in%c("cytolytic", "apm"), 
                                                                                                geometric.mean(meth.data[meth.promoter,x]), 
                                                                                                median(meth.data[meth.promoter,x])),
                                                                                         0)),
                          Meth_Promoter_Var = sapply(colnames(meth.data), function(x) ifelse(length(meth.promoter[which(meth.promoter%in%names(meth.sd))])>0, 
                                                                                             ifelse(length(meth.promoter[which(meth.promoter%in%names(meth.sd))])>1, 
                                                                                                    ifelse(gene.type%in%c("cytolytic",  "apm"),
                                                                                                           geometric.mean(meth.data[meth.promoter[which(meth.promoter%in%names(meth.sd))],x]), 
                                                                                                           median(meth.data[meth.promoter[which(meth.promoter%in%names(meth.sd))],x])), 
                                                                                                    meth.data[meth.promoter[which(meth.promoter%in%names(meth.sd))],x]),
                                                                                             0)),
                          Meth_Body = sapply(colnames(meth.data), function(x) ifelse(length(meth.body)>0, 
                                                                                     ifelse(gene.type%in%c("cytolytic",  "apm"),
                                                                                            geometric.mean(meth.data[meth.body,x]), 
                                                                                            median(meth.data[meth.body,x])),
                                                                                     0)),
                          Meth_Body_Var = sapply(colnames(meth.data), function(x) ifelse(length(meth.body[which(meth.body%in%names(meth.sd))])>0, 
                                                                                         ifelse(length(meth.body[which(meth.body%in%names(meth.sd))])>1, 
                                                                                                ifelse(gene.type%in%c("cytolytic", "apm"),
                                                                                                       geometric.mean(meth.data[meth.body[which(meth.body%in%names(meth.sd))],x]), 
                                                                                                       median(meth.data[meth.body[which(meth.body%in%names(meth.sd))],x])), 
                                                                                                meth.data[meth.body[which(meth.body%in%names(meth.sd))],x]),
                                                                                         0)),
                          RNA_Average = sapply(colnames(meth.data), function(x) ifelse(x%in%colnames(RNA.info), 
                                                                                       ifelse(gene.type%in%c("cytolytic", "apm"),
                                                                                              geometric.mean(log2((RNA.info[,x]+1))), 
                                                                                              median(log2((RNA.info[,x]+1)))),
                                                                                       -10))
                          )
  } else if (!is.null(meth.data)) {
    gene.avg = data.frame(Gene = rep(gene.type, ncol(meth.data)),
                          Sample_Name = colnames(meth.data),
                          Gp_by = sapply(colnames(meth.data),  function(x) samp.info[,col.group][which(samp.info$Sample_Name==x)]),
                          Meth_Promoter = sapply(colnames(meth.data), function(x) ifelse(length(meth.promoter)>0, 
                                                                                         ifelse(gene.type%in%c("cytolytic", "apm"),
                                                                                                geometric.mean(meth.data[meth.promoter,x]), 
                                                                                                median(meth.data[meth.promoter,x])),
                                                                                         0)),
                          Meth_Promoter_Var = sapply(colnames(meth.data), function(x) ifelse(length(meth.promoter[which(meth.promoter%in%names(meth.sd))])>0, 
                                                                                             ifelse(length(meth.promoter[which(meth.promoter%in%names(meth.sd))])>1, 
                                                                                                    ifelse(gene.type%in%c("cytolytic", "apm"),
                                                                                                           geometric.mean(meth.data[meth.promoter[which(meth.promoter%in%names(meth.sd))],x]), 
                                                                                                           median(meth.data[meth.promoter[which(meth.promoter%in%names(meth.sd))],x])), 
                                                                                                    meth.data[meth.promoter[which(meth.promoter%in%names(meth.sd))],x]),
                                                                                             0)),
                          Meth_Body = sapply(colnames(meth.data), function(x) ifelse(length(meth.body)>0, 
                                                                                     ifelse(gene.type%in%c("cytolytic", "apm"),
                                                                                            geometric.mean(meth.data[meth.body,x]), 
                                                                                            median(meth.data[meth.body,x])),
                                                                                     0)),
                          Meth_Body_Var = sapply(colnames(meth.data), function(x) ifelse(length(meth.body[which(meth.body%in%names(meth.sd))])>0, 
                                                                                         ifelse(length(meth.body[which(meth.body%in%names(meth.sd))])>1, 
                                                                                                ifelse(gene.type%in%c("cytolytic", "apm"),
                                                                                                       geometric.mean(meth.data[meth.body[which(meth.body%in%names(meth.sd))],x]), 
                                                                                                       median(meth.data[meth.body[which(meth.body%in%names(meth.sd))],x])), 
                                                                                                meth.data[meth.body[which(meth.body%in%names(meth.sd))],x]),
                                                                                         0))
                          )
  } else if (!is.null(RNA.info)) {
    gene.avg = data.frame(Gene = rep(gene.type, length(grep("_T", colnames(RNA.info)))),
                          Sample_Name = colnames(RNA.info)[grep("_T", colnames(RNA.info))],
                          Gp_by = sapply(colnames(RNA.info)[grep("_T", colnames(RNA.info))],  function(x) samp.info[,col.group][which(samp.info$Sample_Name==x)]),
                          RNA_Average = sapply(colnames(RNA.info)[grep("_T", colnames(RNA.info))], function(x) ifelse(x%in%colnames(RNA.info), 
                                                                                                                      ifelse(gene.type%in%c("cytolytic", "apm"),
                                                                                                                             geometric.mean(log2((RNA.info[,x]+1))), 
                                                                                                                             median(log2((RNA.info[,x]+1)))),
                                                                                                                      -10))
    )
  } else {
    stop("There is no data to plot.")
  }
  if (!is.null(cnv.data)){
    gene.avg %<>% dplyr::mutate(CNV_value = sapply(gene.avg$Sample_Name, function(x) ifelse(length(genes)==1, 
                                                                                           ifelse(gsub("_T", "", x)%in%cnv.data$Sample_Name&genes%in%cnv.data$Genes,
                                                                                                  dplyr::filter(cnv.data, Genes==genes, Sample_Name==gsub("_T", "", x))$CN_values, 
                                                                                                  -10),
                                                                                           -10)))
  } 
  gp.levels = unique(samp.info[,col.group])[order(unique(samp.info[,col.group]))]
  if (!is.null(meth.data)) {
    # create plot of methylation data
    meth.plots = c(ifelse(no.sd.promoter, "Meth_Promoter", "Meth_Promoter_Var"), ifelse(no.sd.body, "Meth_Body", "Meth_Body_Var")) # "Meth_Promoter", 
    meth.grob = list()
    for (meth.plot in meth.plots){
      if (!length(unique(gene.avg[,meth.plot]))==1){
        plot.data = dplyr::mutate(gene.avg, 
                                  Meth_Average=gene.avg[,meth.plot])
        if (col.group=="Survival_Group"){
          plot.data %<>% dplyr::mutate(Treatment = sapply(Sample_Name, function(x) samp.info$Treatment[which(samp.info$Sample_Name==x)]))
        }
        p.meth = ggplot(data=plot.data, aes(x=Gp_by, y=Meth_Average))
        if (violin){
          p.meth = p.meth + geom_violin(width=0.9, aes(fill=Gp_by, color=Gp_by)) + 
            stat_summary(fun=mean, geom="point", size=2, color=cat.cols.bold[[col.group]][gp.levels])
        } else {
          if (col.group=="Survival_Group"){
            p.meth = p.meth + geom_boxplot(outlier.shape=NA, aes(fill=Gp_by, color=Gp_by))
          } else {
            p.meth = p.meth + geom_boxplot(aes(fill=Gp_by, color=Gp_by))
          }
          
        }
        if (col.group=="Survival_Group"){
          p.meth = p.meth + geom_point(position=position_jitterdodge(jitter.width = 1, dodge.width = 0, seed=123), pch=21, aes(fill=Treatment), size=2.5) +
            scale_fill_manual(values = c(cat.cols[[col.group]], cat.cols[["Treatment"]])) + scale_color_manual(values = c(cat.cols.bold[[col.group]], cat.cols.bold[["Treatment"]])) 
        } else {
          p.meth = p.meth + scale_fill_manual(values = cat.cols[[col.group]]) + scale_color_manual(values = cat.cols.bold[[col.group]])
        }
        p.meth = p.meth + 
          labs(title=paste(ifelse(length(genes)==1, genes, paste0(gene.type, " genes")), ifelse(gene.type%in%c("cytolytic", "apm"), "Geometric Mean", "Median"), "Methylation by Group."), 
               y=paste(ifelse(gene.type%in%c("cytolytic", "apm"), "GeoMean", "Median"), ifelse(length(grep("_Var", meth.plot))>0, "most variable", ""),  ifelse(length(grep("_Promoter", meth.plot))>0, "Promoter", "Body"), "Beta Value"), 
               x=gsub("_", " ", gsub("Meth", "Methylation", col.group))) + ylim(0, 1) +
          theme(panel.background=element_rect(fill='transparent'),
                axis.line = element_line(color="black"),
                axis.text.x.bottom = element_text(size=10), 
                axis.text.y.left = element_text(size=10),
                axis.title.x=element_text(size=12),
                axis.title.y=element_text(size=12),
                plot.title = element_text(size=12),
                legend.position = "none") +
          stat_compare_means(aes(group=Gp_by), method=ifelse(length(gp.levels)>2, multi.test, single.test), comparisons=NULL) 
        meth.grob[[meth.plot]] = ggplotGrob(p.meth)
      } else {
        meth.grob[[meth.plot]] = NULL
      }
    }
  } else {
    meth.grob = NULL
  }
  if (!is.null(RNA.info)) {
    # Filter out non rna seq samples
    gene.avg %<>% dplyr::filter(!RNA_Average==-10)
    if (col.group=="Survival_Group"){
      gene.avg %<>% dplyr::mutate(Treatment = sapply(Sample_Name, function(x) samp.info$Treatment[which(samp.info$Sample_Name==x)]))
    }
    # create plot of RNA SEq data
    p.exp = ggplot(data=gene.avg, aes(x=Gp_by, y=RNA_Average))
    if (violin){
      p.exp = p.exp + geom_violin(width=0.9, aes(fill=Gp_by, color=Gp_by)) + 
        stat_summary(fun=mean, geom="point", size=2, color=cat.cols.bold[[col.group]][gp.levels])
    } else {
      if (col.group=="Survival_Group"){
        p.exp = p.exp + geom_boxplot(outlier.shape = NA,  aes(fill=Gp_by, color=Gp_by))
      } else {
        p.exp = p.exp + geom_boxplot(aes(fill=Gp_by, color=Gp_by))
      }
    }
    if (col.group=="Survival_Group"){
      p.exp = p.exp + geom_point(position=position_jitterdodge(jitter.width = 1, dodge.width = 0, seed=123), pch=21, aes(fill=Treatment), size = 2.5) +
        scale_fill_manual(values = c(cat.cols[[col.group]], cat.cols[["Treatment"]])) + scale_color_manual(values = c(cat.cols.bold[[col.group]], cat.cols.bold[["Treatment"]])) 
    } else {
      p.exp = p.exp + scale_fill_manual(values = cat.cols[[col.group]]) + scale_color_manual(values = cat.cols.bold[[col.group]])
    }
    p.exp = p.exp + 
      labs(title="",
           y=ifelse(gene.type%in%c("cytolytic", "apm"), "APM \nGeometric Mean Log2(exp+1)", ifelse(length(genes)==1, paste(genes, "Log2(exp+1)"), "Median Log2(exp+1)")), 
           x="") + 
      theme(panel.background=element_rect(fill='transparent'),
            axis.line = element_line(color="black"),
            axis.text.x.bottom = element_text(size=14), 
            axis.text.y.left = element_text(size=12),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=16),
            plot.title = element_text(size=12),
            legend.position = "none")
    if (all.comp){
      my.comparisons = list(c("Group1", "Group2"), c("Group1", "Group3"), c("Group1", "Group4"), c("Group2", "Group3"), c("Group2", "Group4"), c("Group3", "Group4"))
      p.exp = p.exp + stat_compare_means(method=single.test, comparisons = my.comparisons, label.y = c(max(gene.avg$RNA_Average)+max(gene.avg$RNA_Average)*0.08, 
                                                                                                       max(gene.avg$RNA_Average)+max(gene.avg$RNA_Average)*0.16, 
                                                                                                       max(gene.avg$RNA_Average)+max(gene.avg$RNA_Average)*0.24,
                                                                                                       max(gene.avg$RNA_Average)+max(gene.avg$RNA_Average)*0.32,
                                                                                                       max(gene.avg$RNA_Average)+max(gene.avg$RNA_Average)*0.40,
                                                                                                       max(gene.avg$RNA_Average)+max(gene.avg$RNA_Average)*0.48)) + 
        stat_compare_means(aes(group=Gp_by), method=ifelse(length(gp.levels)>2, multi.test, single.test), comparisons=NULL) + 
        ylim(0, max(gene.avg$RNA_Average)+max(gene.avg$RNA_Average)*0.55)
    } else {
      p.exp = p.exp + stat_compare_means(aes(group=Gp_by), method=ifelse(length(gp.levels)>2, multi.test, single.test), comparisons=NULL) 
    }
    exp.grob = ggplotGrob(p.exp)
  } else {
    exp.grob = NULL
  }
  if (!is.null(meth.data) & !is.null(RNA.info)) {
    scat.grob = list()
    scat.plots = c("Meth_Promoter", "Meth_Promoter_Var", "Meth_Body", "Meth_Body_Var", "CNV_value")
    # create correlation plot between Methylation and RNASeq
    for (scat.plot in scat.plots){
      if (scat.plot %in% colnames(gene.avg)){
        scat = ggscatter(gene.avg, 
                         x=scat.plot, 
                         y="RNA_Average", 
                         color = "Gp_by",
                         fill= "Gp_by",
                         size=2,
                         add="reg.line", 
                         add.params=list(color="darkgrey", fill="#eeeeee"),
                         conf.int=TRUE,
                         cor.coef=TRUE,
                         cor.coef.coord = c((min(gene.avg[,scat.plot])+max(gene.avg[,scat.plot]))/2, max(gene.avg$RNA_Average)+1),
                         cor.method="pearson", 
                         title=paste(ifelse(length(genes)==1, genes, paste0(gene.type, " genes")), "Expression vs ", ifelse(length(grep("Var", scat.plot))>0, paste0("most variable ", gsub("Meth_", "", gsub("_Var", "", scat.plot)), " Beta Values."), ifelse(!scat.plot=="CNV_value", paste0(gsub("Meth_", " ", scat.plot), " Beta Values."), "Raw Copy Number"))),
                         xlab=paste(ifelse(gene.type%in%c("cytolytic", "apm"), "GeoMean Beta Value", ifelse(scat.plot=="CNV_value", "Raw Copy Number", "Median Beta Value"))),
                         ylab=paste(ifelse(gene.type%in%c("cytolytic", "apm"), "GeoMean", "Median"), "Log2(exp+1)"),
                         ylim=c(0,max(gene.avg$RNA_Average)+1), # min(gene.avg$RNA_Average)
                         xlim=c(0, max(gene.avg[,scat.plot])), # min(gene.avg[,scat.plot])
                         palette=cat.cols.bold[[col.group]]) 
        scat = scat + theme(plot.title = element_text(size = 12), 
                            axis.title.x=element_text(size=12),
                            axis.title.y=element_text(size=12),
                            axis.text.x.bottom = element_text(size=10),
                            axis.text.y.left = element_text(size=10),
                            legend.position = "none")
        scat.grob[[scat.plot]] = ggplotGrob(scat)
      } else {
        scat.grob[[scat.plot]] = NULL
      }
    }
  } else {
    scat.grob = NULL
  }
  plot.list = list(data=gene.avg, meth=meth.grob, rna=exp.grob, corr=scat.grob)
  return(plot.list)
}
