#' Calculate data attribute statistics
#' @description
#' Calculates the similarity of data attributes between the first matrix and the other matrix.
#'
#'
#' @param maxNForDisp The maximal number of samples that will be used to estimate dispersions. By default, all samples are used. This can be lowered to speed up calculations (and obtain approximate results) for large data sets.
#' @param maxNForCorr The maximal number of samples (features) for which pairwise correlation coefficients will be calculated. If the number of samples (features) exceeds this number, they will be randomly subsampled.
#' @param subsampleSize The number of randomly selected observations (samples, features or pairs of samples or features) for which certain (time-consuming) statistics will be calculated.
#' @param kfrac For statistics that require the extraction of the k nearest neighbors of a given point, the number of neighbors will be max(kmin, kfrac * nrow(df))
#' @param kmin For statistics that require the extraction of the k nearest neighbors of a given point, the number of neighbors will be max(kmin, kfrac * nrow(df))
#' @param matList List of gene expression count matrices
#'
#' @return Calculate data attribute statistics
#'
#'
#' @export
calcDataProperty <- function(matList, maxNForDisp, maxNForCorr,
                  subsampleSize, kfrac, kmin) {
  set.seed(123)
  print("Calculate data attribute statistics")

  if (!is(matList, "list")) {
    stop("matList must be a list.", call. = FALSE)
  }
  if (length(setdiff(unique(names(matList)),
                     c("", NA, NULL))) != length(matList)) {
    stop("matList must be a named list, ",
         "with a unique name for each element.", call. = FALSE)
  }
  if (!all(vapply(matList, function(w) {
    is(w, "data.frame") | is(w, "matrix")
  }, FALSE))) {
    stop("All elements of matList must be data.frames ",
         "or matrices. See the DESeq2 Bioconductor package ",
         "(http://bioconductor.org/packages/release/bioc/html/DESeq2.html) ",
         "for more information about the DESeqDataSet class.", call. = FALSE)
  }

  ## If some objects are data frames or matrices, convert them into
  ## DESeqDataSets
  matList <- lapply(matList, function(ds) {
    if (is(ds, "DESeqDataSet")) {
      ds
    } else {
      DESeq2::DESeqDataSetFromMatrix(
        countData = round(as.matrix(ds)),
        colData = data.frame(sample = seq_len(ncol(ds))),
        design = ~ 1)
    }
  })
  stopifnot(all(vapply(matList, function(w) is(w, "DESeqDataSet"), FALSE)))



  ##数据集数
  nDatasets <- length(matList)
  message(paste("There are a total of ", nDatasets, "datasets"))

  if(nDatasets==1){
    stop("matList must be a named list, ",
         "must have at least two elements.", call. = FALSE)
  }

  message("Calculating dispersions. This can take some time, please be patient.")

  if (is.finite(maxNForDisp)) {
    msg <- paste0("If there are more than ", maxNForDisp,
                  " samples in a data set, the dispersions are calculated using ",
                  maxNForDisp, " randomly selected samples.")
  } else {
    msg <- ""
  }
  obj <- calculateDispersionsddsList(ddsList = matList,
                                     maxNForDisp = maxNForDisp)
  message("Calculating sample correlations.")
  sampleCorrDF <- calculateSampleCorrs(ddsList = obj,
                                       maxNForCorr = maxNForCorr)
  message("Calculating feature correlations.")
  featureCorrDF <- calculateFeatureCorrs(ddsList = obj,
                                         maxNForCorr = maxNForCorr)

  #细胞属性值计算
  sampleDF <- lapply(obj, function(x) {
    data.frame(
      Libsize = colSums(x$dge$counts),
      Fraczero = colMeans(x$dge$counts == 0),
      TMM = x$dge$samples$norm.factors
    ) %>% dplyr::mutate(EffLibsize = Libsize * TMM)
  })
  ns <- sapply(sampleDF, nrow)
  sampleDF <- do.call(rbind, sampleDF) %>%
    dplyr::mutate(dataset = rep(names(sampleDF), ns))



  #基因属性值计算
  featureDF <- lapply(obj, function(x) {
    data.frame(
      Tagwise = sqrt(x$dge$tagwise.dispersion),
      Common = sqrt(x$dge$common.dispersion),
      Trend = sqrt(x$dge$trended.dispersion),
      AveLogCPM = x$dge$AveLogCPM,
      AveLogCPMDisp = x$dge$AveLogCPMDisp,
      average_log2_cpm = apply(edgeR::cpm(x$dge, prior.count = 2, log = TRUE), 1, mean),
      variance_log2_cpm = apply(edgeR::cpm(x$dge, prior.count = 2, log = TRUE), 1, var),
      Fraczero = rowMeans(x$dge$counts == 0),
      dispGeneEst = rowData(x$dds)$dispGeneEst,
      dispFit = rowData(x$dds)$dispFit,
      dispFinal = rowData(x$dds)$dispersion,
      baseMeanDisp = rowData(x$dds)$baseMeanDisp,
      baseMean = rowData(x$dds)$baseMean
    )
  })
  ns <- sapply(featureDF, nrow)
  featureDF <- do.call(rbind, featureDF) %>%
    dplyr::mutate(dataset = rep(names(featureDF), ns))






  datasetDF <- do.call(rbind, lapply(obj, function(x) {
    data.frame(
      prior_df = paste0("prior.df = ", round(x$dge$prior.df, 2)),
      nVars = nrow(x$dge$counts),
      nSamples = ncol(x$dge$counts)
    )
  })) %>%
    dplyr::mutate(AveLogCPMDisp = 0.8 * max(featureDF$AveLogCPMDisp)) %>%
    dplyr::mutate(Tagwise = 0.9 * max(featureDF$Tagwise)) %>%
    dplyr::mutate(dataset = names(obj))


  message("Calculating data property.")
  AveLogCPM_Tagwise<- makeDF(df = featureDF, column = c("AveLogCPMDisp", "Tagwise"),
                             subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)


  Mean_dispersion <- makeDF(df = featureDF, column = c("baseMeanDisp", "dispGeneEst"),
                            subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)

  Mean_Var <- makeDF(df = featureDF, column = c("average_log2_cpm", "variance_log2_cpm"),
                     subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)

  Lib <- makeDF(df = sampleDF, column = "Libsize",
                subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)
  TMM <- makeDF(df = sampleDF, column = "TMM",
                subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)
  EffLibsize <- makeDF(df = sampleDF, column = "EffLibsize",
                       subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)
  AveLogCPM <- makeDF(df = featureDF, column = "AveLogCPM",
                      subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)
  sampleFraczero <- makeDF(df = sampleDF, column = "Fraczero",
                           subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)
  featureFraczero <- makeDF(df = featureDF, column = "Fraczero",
                            subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)
  sampleCorrelation <- makeDF(df = sampleCorrDF, column = "Correlation",
                              subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)
  featureCorrelation <- makeDF(df = featureCorrDF, column = "Correlation",
                               subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)
  Lib_Fraczero <- makeDF(df = sampleDF, column = c("Libsize", "Fraczero"),
                         subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)
  AveLogCPM_Fraczero <- makeDF(df = featureDF, column = c("AveLogCPM", "Fraczero"),
                               subsampleSize = subsampleSize, kmin = kmin, kfrac = kfrac)



  final_list <- list(AveLogCPM_Tagwise=AveLogCPM_Tagwise, Mean_dispersion=Mean_dispersion,
                     Mean_Var = Mean_Var, Lib=Lib, TMM=TMM, EffLibsize=EffLibsize,
                     AveLogCPM=AveLogCPM, sampleFraczero=sampleFraczero,featureFraczero=featureFraczero,
                     sampleCorrelation=sampleCorrelation,featureCorrelation=featureCorrelation,
                     Lib_Fraczero=Lib_Fraczero,AveLogCPM_Fraczero=AveLogCPM_Fraczero)

  #write.csv(final_list[1], "final.csv", row.names = FALSE)

  return(final_list)
}



StandardizedCalculationScore <- function(tmp_list){

  #standardization
  Runs_statistic <- lapply(tmp_list, function(x) {

    if(!is.null(x$`Runs statistic`)){

      if(all(x$`Runs statistic` == x$`Runs statistic`[1])){
        x$`Runs statistic`=tibble(V1=rep(0, length(x$`Runs statistic`)))
      } else{

        x$`Runs statistic` <- scale(x$`Runs statistic`)  # 使用scale函数进行标准化

        min_val <- min(x$`Runs statistic`,na.rm = TRUE)
        max_val <- max(x$`Runs statistic`,na.rm = TRUE)
        normalized_data <- (x$`Runs statistic` - min_val) / (max_val - min_val)

        x$`Runs statistic` <- as.data.frame(1 - normalized_data)
      }
    }
  })


  Multi_KS <- lapply(tmp_list, function(x) {
    if(!is.null(x$`Multivariate KS statistic`)){

      if(all(x$`Multivariate KS statistic` == x$`Multivariate KS statistic`[1])){
        x$`Multivariate KS statistic`=tibble(V1=rep(0, length(x$`Multivariate KS statistic`)))
      }
      else{

        x$`Multivariate KS statistic` <- abs(x$`Multivariate KS statistic`)   #多元KS需计算绝对值
        x$`Multivariate KS statistic` <- scale(x$`Multivariate KS statistic`)  # 使用scale函数进行标准化

        min_val <- min(x$`Multivariate KS statistic`,na.rm = TRUE)
        max_val <- max(x$`Multivariate KS statistic`,na.rm = TRUE)
        normalized_data <- (x$`Multivariate KS statistic` - min_val) / (max_val - min_val)

        x$`Multivariate KS statistic` <- as.data.frame(1 - normalized_data)
      }

    }



  })
  KDE_T <- lapply(tmp_list, function(x) {
    if(!is.null(x$`KDE T-statistic`)){

      if(all(x$`KDE T-statistic` == x$`KDE T-statistic`[1])){
        x$`KDE T-statistic`=tibble(V1=rep(0, length(x$`KDE T-statistic`)))
      }
      else{
        x$`KDE T-statistic` <- abs(x$`KDE T-statistic`)
        x$`KDE T-statistic` <- scale(x$`KDE T-statistic`)  # 使用scale函数进行标准化

        min_val <- min(x$`KDE T-statistic`,na.rm = TRUE)
        max_val <- max(x$`KDE T-statistic`,na.rm = TRUE)
        normalized_data <- (x$`KDE T-statistic` - min_val) / (max_val - min_val)

        x$`KDE T-statistic` <- as.data.frame(1 - normalized_data)
      }

    }


  })

  KDE_Z <- lapply(tmp_list, function(x) {


    if(!is.null(x$`KDE Z-statistic`)){

      if(all(x$`KDE Z-statistic` == x$`KDE Z-statistic`[1])){
        x$`KDE Z-statistic`=tibble(V1=rep(0, length(x$`KDE Z-statistic`)))
      }
      else{
        x$`KDE Z-statistic` <- abs(x$`KDE Z-statistic`)
        x$`KDE Z-statistic` <- scale(x$`KDE Z-statistic`)  # 使用scale函数进行标准化

        min_val <- min(x$`KDE Z-statistic`,na.rm = TRUE)
        max_val <- max(x$`KDE Z-statistic`,na.rm = TRUE)
        normalized_data <- (x$`KDE Z-statistic` - min_val) / (max_val - min_val)

        x$`KDE Z-statistic` <- as.data.frame(1 - normalized_data)
      }
    }
  })


  tmp_list$AveLogCPM_Tagwise$`Multivariate KS statistic` <- Multi_KS$AveLogCPM_Tagwise$V1
  tmp_list$AveLogCPM_Tagwise$`KDE T-statistic` <- KDE_T$AveLogCPM_Tagwise$V1
  tmp_list$AveLogCPM_Tagwise$`KDE Z-statistic` <- KDE_Z$AveLogCPM_Tagwise$V1



  tmp_list$Mean_dispersion$`Multivariate KS statistic` <- Multi_KS$Mean_dispersion$V1
  tmp_list$Mean_dispersion$`KDE T-statistic` <- KDE_T$Mean_dispersion$V1
  tmp_list$Mean_dispersion$`KDE Z-statistic` <- KDE_Z$Mean_dispersion$V1


  tmp_list$Mean_Var$`Multivariate KS statistic` <- Multi_KS$Mean_Var$V1
  tmp_list$Mean_Var$`KDE T-statistic` <- KDE_T$Mean_Var$V1
  tmp_list$Mean_Var$`KDE Z-statistic` <- KDE_Z$Mean_Var$V1


  tmp_list$Lib_Fraczero$`Multivariate KS statistic` <- Multi_KS$Lib_Fraczero$V1
  tmp_list$Lib_Fraczero$`KDE T-statistic` <- KDE_T$Lib_Fraczero$V1
  tmp_list$Lib_Fraczero$`KDE Z-statistic` <- KDE_Z$Lib_Fraczero$V1

  tmp_list$AveLogCPM_Fraczero$`Multivariate KS statistic` <- Multi_KS$AveLogCPM_Fraczero$V1
  tmp_list$AveLogCPM_Fraczero$`KDE T-statistic` <- KDE_T$AveLogCPM_Fraczero$V1
  tmp_list$AveLogCPM_Fraczero$`KDE Z-statistic` <- KDE_Z$AveLogCPM_Fraczero$V1

  tmp_list$Lib$`Runs statistic` <- Runs_statistic$Lib$V1
  tmp_list$TMM$`Runs statistic` <- Runs_statistic$TMM$V1
  tmp_list$EffLibsize$`Runs statistic` <- Runs_statistic$EffLibsize$V1
  tmp_list$AveLogCPM$`Runs statistic` <- Runs_statistic$AveLogCPM$V1
  tmp_list$sampleFraczero$`Runs statistic` <- Runs_statistic$sampleFraczero$V1
  tmp_list$featureFraczero$`Runs statistic` <- Runs_statistic$featureFraczero$V1
  tmp_list$sampleCorrelation$`Runs statistic` <- Runs_statistic$sampleCorrelation$V1
  tmp_list$featureCorrelation$`Runs statistic` <- Runs_statistic$featureCorrelation$V1

  #Calculate Score
  tmp_list <- lapply(tmp_list, function(element) {

    if (length(colnames(element)) == 8) {
      element$score <- (1 - element$`NN rejection fraction`) +
        (1 - abs(element$`Average silhouette width`))/2 +
        (1 - abs(element$`Average local silhouette width`))/2 +
        element$`Multivariate KS statistic` +
        element$`KDE Z-statistic`
    }else{
      element$score <- (1-element$`K-S statistic`)+(1-element$`Scaled area between eCDFs`)+
        (1-element$`Runs statistic`)+(1-element$`NN rejection fraction`)+
        (1-abs(element$`Average silhouette width`))/2+(1-abs(element$`Average local silhouette width`))/2
    }
    return(element)
  })
  rank <- tmp_list$AveLogCPM_Tagwise[c('dataset1','dataset2')]
  rank$score <- 0


  #sum Score
  for(element in tmp_list){
    for(i in 1:nrow(rank)){
      rank$score[i]<-element$score[i]+rank$score[i]
    }
  }

  #rank sum Score
  sorted_rank <- rank[order(-rank$score), ]
  rownames(sorted_rank) <- rownames(rank)

  #rank the scores of each attribute
  tmp_list <- lapply(tmp_list, function(element) {
    sorted_element <- element[order(element$score, decreasing = TRUE), ]
  })

  tmp_list$sorted_rank <- sorted_rank
  return(tmp_list)

}



#' Calculate data attribute statistics and normalized data
#' @description
#' This method calculates the similarity of data attributes between the first SCE object and the other SCE objects, standardizes the similarity scores, and obtains a similarity ranking.
#'
#' @param SCEList List of SingleCellExperiment
#' @param maxNForDisp The maximal number of samples that will be used to estimate dispersions. By default, all samples are used. This can be lowered to speed up calculations (and obtain approximate results) for large data sets.
#' @param maxNForCorr The maximal number of samples (features) for which pairwise correlation coefficients will be calculated. If the number of samples (features) exceeds this number, they will be randomly subsampled.
#' @param subsampleSize The number of randomly selected observations (samples, features or pairs of samples or features) for which certain (time-consuming) statistics will be calculated.
#' @param kfrac For statistics that require the extraction of the k nearest neighbors of a given point, the number of neighbors will be max(kmin, kfrac * nrow(df))
#' @param kmin For statistics that require the extraction of the k nearest neighbors of a given point, the number of neighbors will be max(kmin, kfrac * nrow(df))
#'
#' @import SingleCellExperiment
#'
#' @return Calculate data attribute statistics and normalized data
#'
#' @examples
#' data("ESCO_sce_293T_jurkat", "muscat_sce_293T_jurkat", "scDesign3_sce_293T_jurkat", "origin_sce_293T_jurkat")
#' SCEList <- list(origin=origin_sce_293T_jurkat, ESCO=ESCO_sce_293T_jurkat, muscat=muscat_sce_293T_jurkat, scDesign3 = scDesign3_sce_293T_jurkat)
#' test <- calcDataPropertyScores(SCEList)
#'
#' @export
calcDataPropertyScores <- function(SCEList, maxNForDisp = Inf, maxNForCorr = 25,
                                   subsampleSize = 1000, kfrac = 0.05, kmin = 5){

  counts_list <- lapply(SCEList, function(sce) {
    # 提取 counts 矩阵
    counts <- assays(sce)$counts
    return(counts)
  })

  # 保留原始名称（如果 SCEList 是命名列表）
  names(counts_list) <- names(SCEList)

  result_list <- calcDataProperty(matList = counts_list, maxNForDisp = maxNForDisp, maxNForCorr = maxNForCorr,
                                  subsampleSize = subsampleSize, kfrac = kfrac, kmin = kmin)

  normalized_list <- StandardizedCalculationScore(result_list)

  fin_list <- list(result_list = result_list, normalized_list = normalized_list)

  return(fin_list)

}

