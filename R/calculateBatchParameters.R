#' Calculate clustering evaluation indicators for simulated data
#'
#' @param mat Gene expression matrix
#' @param batch_info Cell Batch Label
#' @param TSNE_ntop A numerical scalar that specifies the number of features with the highest variance used for TSNE dimensionality reduction.
#' @param PCA_ntop A numerical scalar that specifies the number of features with the highest variance used for PCA dimensionality reduction.
#'
#' @return The value or list of clustering evaluation
#'
#'
#' @export
#'
calculateBatch <- function(mat, batch_info, TSNE_ntop, PCA_ntop){

  # import mclust stats scater kBET

  #将ddsList提取并赋值为matrix
  matrix <- as.matrix(mat)

  if (ncol(matrix) != length(batch_info)) {
    stop("Error: The number of columns in the matrix does not match the length of batch_info.")
  }

  #转换为SCE对象
  sce <- SingleCellExperiment(
    assays = list(counts = matrix)
  )

  sce <- scater::logNormCounts(sce)

  # 执行降维
  message("reducing the dimensionality of the raw data.")
  sce <- scater::runTSNE(sce, ntop = TSNE_ntop)
  sce <- scater::runPCA(sce, ntop = PCA_ntop)
  #sce <- runUMAP(sce)

  tsne_matrix <- SingleCellExperiment::reducedDim(sce, "TSNE")
  pca_matrix <- SingleCellExperiment::reducedDim(sce, "PCA")
  #umap_matrix <- reducedDim(sce, "UMAP")

  tsne_dist <- stats::dist(tsne_matrix)
  pca_dist <- stats::dist(pca_matrix)
  #umap_dist <- stats::dist(umap_matrix)

  k <- round(sqrt(length(batch_info)))

  #kBET计算
  message("Calculating kBET.")
  pca_kBET <- kBET::kBET(df = t(matrix), k0 = k, batch = as.numeric(as.factor(batch_info)), do.pca = TRUE, plot = FALSE)
  tsne_kBET <- kBET::kBET(df = tsne_matrix, k0 = k, batch = as.numeric(as.factor(batch_info)), do.pca = FALSE, plot = FALSE)
  #umap_kBET <- kBET(df = umap_matrix, k0 = k, batch = as.numeric(as.factor(batch_info)), do.pca = FALSE, plot = FALSE)

  #batch silhouette计算
  message("Calculating ASW.")
  pca_silhouette_width <- cluster::silhouette(x = as.numeric(as.factor(batch_info)),
                                          pca_dist)
  pca_average_silhouette <- mean(pca_silhouette_width[, 3])

  tsne_silhouette_width <- cluster::silhouette(x = as.numeric(as.factor(batch_info)),
                                               tsne_dist)
  tsne_average_silhouette <- mean(tsne_silhouette_width[, 3])

  #umap_silhouette_width <- cluster::silhouette(x = as.numeric(as.factor(batch_info)),
  #                                             umap_dist)
  #umap_average_silhouette <- mean(umap_silhouette_width[, 3])

  #LISI
  message("Calculating LISI.")
  pca_LISI <- lisi::compute_lisi(pca_matrix,
                                 meta_data =  data.frame("batch" = batch_info),
                                 label_colnames = "batch",
                                 perplexity = k)
  pca_meanLISI <- mean(pca_LISI$batch)

  tsne_LISI <- lisi::compute_lisi(tsne_matrix,
                                 meta_data =  data.frame("batch" = batch_info),
                                 label_colnames = "batch",
                                 perplexity = k)
  tsne_meanLISI <- mean(tsne_LISI$batch)

  #umap_LISI <- lisi::compute_lisi(umap_matrix,
  #                               meta_data =  data.frame("batch" = batch_info),
  #                               label_colnames = "batch",
  #                               perplexity = k)
  #umap_meanLISI <- mean(umap_LISI$batch)


  #ARI
  message("Calculating ARI.")
  pca_kmeans_result <- kmeans(pca_matrix, centers = 2, nstart = 25)
  pca_cluster_labels <- pca_kmeans_result$cluster
  pca_ARI <- mclust::adjustedRandIndex(pca_cluster_labels,batch_info)

  tsne_kmeans_result <- kmeans(tsne_matrix, centers = 2, nstart = 25)
  tsne_cluster_labels <- tsne_kmeans_result$cluster
  tsne_ARI <- mclust::adjustedRandIndex(tsne_cluster_labels,batch_info)

  #umap_kmeans_result <- kmeans(umap_matrix, centers = 2, nstart = 25)
  #umap_cluster_labels <- umap_kmeans_result$cluster
  #umap_ARI <- mclust::adjustedRandIndex(umap_cluster_labels,batch_info)

  data_list <- list(
    kBET = list(
      PCA_kBET = pca_kBET,
      TSNE_kBET = tsne_kBET
    ),
    ASW = list(
      PCA_ASW = pca_average_silhouette,
      TSNE_ASW = tsne_average_silhouette
    ),
    LISI = list(
      PCA_LISI = pca_meanLISI,
      TSNE_LISI = tsne_meanLISI
    ),
    ARI = list(
      PCA_ARI = pca_ARI,
      TSNE_ARI = tsne_ARI
    )
  )
  return(data_list)
}

#' Calculate batch effect evaluation indicators and normalize.
#'
#' @description
#' This method calculates batch effects for each SCE object and normalizes the indicators. The higher the normalization value, the more prominent the batch effect.
#'
#' @param SCEList SCEList is a list containing multiple SingleCellExperiment objects.Please ensure that each SCE object has a 'Batch' column in its colData.
#' @param TSNE_ntop A numerical scalar that specifies the number of features with the highest variance used for TSNE dimensionality reduction.
#' @param PCA_ntop A numerical scalar that specifies the number of features with the highest variance used for PCA dimensionality reduction.
#'
#' @examples
#' data("splatter_sce_Batch", "SPARSim_sce_Batch")
#' SCEList <- list(SPARSim = SPARSim_sce_Batch, splatter=splatter_sce_Batch)
#' test <- calculateBatchParameters(SCEList)
#'
#'
#' @return Batch effect evaluation list and normalized values
#' @export

calculateBatchParameters <- function(SCEList, TSNE_ntop = 500, PCA_ntop = 500){

  if (!is(SCEList, "list")) {
    stop("SCEList must be a list.", call. = FALSE)
  }

  if (length(setdiff(unique(names(SCEList)),
                     c("", NA, NULL))) != length(SCEList)) {
    stop("If SCEList is List, it must be a named list, ",
         "with a unique name for each element.", call. = FALSE)
  }

  ##数据集数
  nDatasets <- length(SCEList)
  message(paste("There are a total of ", nDatasets, "datasets"))


  result_list <- lapply(SCEList, function(sce){
    expr_data <- assays(sce)$counts
    Batch <- sce$Batch
    calculateBatch(expr_data, Batch, TSNE_ntop = TSNE_ntop, PCA_ntop = PCA_ntop)
  })

  tmp_list <- vector("list", length(result_list))
  names(tmp_list) <- names(result_list)  # 保留原名称

  for (i in seq_along(result_list)) {
    # 提取当前元素的值
    current_pca_kBET <- result_list[[i]]$kBET$PCA_kBET$summary$kBET.observed[1]
    current_tsne_kBET <- result_list[[i]]$kBET$TSNE_kBET$summary$kBET.observed[1]

    current_pca_asw <- result_list[[i]]$ASW$PCA_ASW
    current_tsne_asw <- result_list[[i]]$ASW$TSNE_ASW

    current_pca_LISI <- result_list[[i]]$LISI$PCA_LISI
    current_tsne_LISI <- result_list[[i]]$LISI$TSNE_LISI

    current_pca_ARI <- result_list[[i]]$ARI$PCA_ARI
    current_tsne_ARI <- result_list[[i]]$ARI$TSNE_ARI

    tmp_list[[i]] <- list(
      kBET = list(
        PCA_kBET = current_pca_kBET,
        TSNE_kBET = current_tsne_kBET
        #UMAP_asw = diff_umap_asw
      ),
      ASW = list(
        PCA_ASW = current_pca_asw,
        TSNE_ASW = current_tsne_asw
      ),
      LISI = list(
        PCA_LISI = current_pca_LISI,
        TSNE_LISI = current_tsne_LISI
      ),
      ARI = list(
        PCA_ARI = current_pca_ARI,
        TSNE_ARI = current_tsne_ARI
      )
    )
  }

  # 标准化与归一化函数
  scale_and_normalize <- function(data_matrix, center = TRUE, scale = TRUE) {

    if(all(data_matrix == data_matrix[1])){
      normalized_data=array(1, dim = length(data_matrix))
    }
    else{
      # 1. 标准化数据
      scaled_data <- scale(data_matrix, center = center, scale = scale)

      # 2. 计算最小最大值（忽略NA）
      min_val <- min(scaled_data, na.rm = TRUE)
      max_val <- max(scaled_data, na.rm = TRUE)

      # 3. 最小-最大归一化到[0,1]范围
      normalized_data <- (scaled_data - min_val) / (max_val - min_val)
    }


    # 保留原始的行列名
    rownames(normalized_data) <- names(data_matrix)

    return(normalized_data)
  }

  #获得对应指标的值，并将其归一化
  pca_kBET_matrix <- sapply(tmp_list, function(x) abs(x$kBET$PCA_kBET))
  pca_kBET_normalize <- scale_and_normalize(pca_kBET_matrix)

  tsne_kBET_matrix <- sapply(tmp_list, function(x) abs(x$kBET$TSNE_kBET))
  tsne_kBET_normalize <- scale_and_normalize(tsne_kBET_matrix)

  pca_ASW_matrix <- sapply(tmp_list, function(x) abs(x$ASW$PCA_ASW))
  pca_ASW_normalize <- scale_and_normalize(pca_ASW_matrix)

  tsne_ASW_matrix <- sapply(tmp_list, function(x) abs(x$ASW$TSNE_ASW))
  tsne_ASW_normalize <- scale_and_normalize(tsne_ASW_matrix)

  pca_LISI_matrix <- sapply(tmp_list, function(x) -abs(x$LISI$PCA_LISI))
  pca_LISI_normalize <- scale_and_normalize(pca_LISI_matrix)

  tsne_LISI_matrix <- sapply(tmp_list, function(x) -abs(x$LISI$TSNE_LISI))
  tsne_LISI_normalize <- scale_and_normalize(tsne_LISI_matrix)

  pca_ARI_matrix <- sapply(tmp_list, function(x) abs(x$ARI$PCA_ARI))
  pca_ARI_normalize <- scale_and_normalize(pca_ARI_matrix)

  tsne_ARI_matrix <- sapply(tmp_list, function(x) abs(x$ARI$TSNE_ARI))
  tsne_ARI_normalize <- scale_and_normalize(tsne_ARI_matrix)

  score_list <- tmp_list

  #将数据放入score_list
  for (i in seq_along(tmp_list)) {
    score_list[[i]]$kBET$PCA_kBET <- pca_kBET_normalize[[i]]
    score_list[[i]]$kBET$TSNE_kBET <- tsne_kBET_normalize[[i]]

    score_list[[i]]$ASW$PCA_ASW <- pca_ASW_normalize[[i]]
    score_list[[i]]$ASW$TSNE_ASW <- tsne_ASW_normalize[[i]]

    score_list[[i]]$LISI$PCA_LISI <- pca_LISI_normalize[[i]]
    score_list[[i]]$LISI$TSNE_LISI <- tsne_LISI_normalize[[i]]

    score_list[[i]]$ARI$PCA_ARI <- pca_ARI_normalize[[i]]
    score_list[[i]]$ARI$TSNE_ARI <- tsne_ARI_normalize[[i]]
  }
  #计算PCA和TSNE和sum
  for (i in seq_along(score_list)) {
    score_list[[i]]$score$PCA <- (score_list[[i]]$ASW$PCA_ASW + score_list[[i]]$kBET$PCA_kBET + score_list[[i]]$LISI$PCA_LISI + score_list[[i]]$ARI$PCA_ARI)/4

    score_list[[i]]$score$TSNE <- (score_list[[i]]$ASW$TSNE_ASW + score_list[[i]]$kBET$TSNE_kBET + score_list[[i]]$LISI$TSNE_LISI + score_list[[i]]$ARI$TSNE_ARI)/4

    score_list[[i]]$score$sum <- (score_list[[i]]$score$PCA + score_list[[i]]$score$TSNE)/2
  }

  sum <- sapply(score_list, function(x) x$score$sum)

  result_df <- data.frame(
    sum_score = sum,
    row.names = names(score_list)  # 设置行名为列表名
  )

  score_list$score <- result_df

  return(list(result_list = result_list, normalized_list = score_list))

}
