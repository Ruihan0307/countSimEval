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

#' Calculate batch effects and return the results as a list
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
#' @return List of Batch Effect evaluation
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
}
