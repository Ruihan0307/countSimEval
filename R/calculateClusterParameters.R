#' Calculate clustering evaluation indicators
#'
#' @param cluster_info Cell group label
#' @param TSNE_ntop A numerical scalar that specifies the number of features with the highest variance used for TSNE dimensionality reduction.
#' @param PCA_ntop A numerical scalar that specifies the number of features with the highest variance used for PCA dimensionality reduction.
#' @param mat Gene expression matrix
#'
#' @return The value or list of clustering evaluation
#'
#'
#' @export
calculateCluster <- function(mat, cluster_info, TSNE_ntop, PCA_ntop){

    #将ddsList提取并赋值为matrix
    matrix <- as.matrix(mat)

    if (ncol(matrix) != length(cluster_info)) {
      stop("Error: The number of columns in the matrix does not match the length of cluster_info.")
    }

    #转换为SCE对象
    sce <- SingleCellExperiment(
      assays = list(counts = matrix)
    )

    sce <- scater::logNormCounts(sce)

    ##使用原始数据计算距离
    #message("使用原始数据计算距离.")
    #raw_dist <- parallelDist::parDist(t(assays(sce)$counts))

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

    message("calculating average silhouette width.")
    #计算ASW
    tsne_silhouette_width <- cluster::silhouette(x = as.numeric(as.factor(cluster_info)),
                                                 tsne_dist)
    tsne_average_silhouette <- mean(tsne_silhouette_width[, 3])

    pca_silhouette_width <- cluster::silhouette(x = as.numeric(as.factor(cluster_info)),
                                                pca_dist)
    pca_average_silhouette <- mean(pca_silhouette_width[, 3])

    #umap_silhouette_width <- cluster::silhouette(x = as.numeric(as.factor(cluster_info)),
    #                                             umap_dist)
    #umap_average_silhouette <- mean(umap_silhouette_width[, 3])

    #raw_silhouette_width <- cluster::silhouette(x = as.numeric(as.factor(cluster_info)),
    #                                            raw_dist)
    #raw_average_silhouette <- mean(raw_silhouette_width[, 3])

    message("calculating Connectivity.")
    #计算连通性
    tsne_con <- clValid::connectivity(distance = tsne_dist, clusters = cluster_info)
    pca_con <- clValid::connectivity(distance = pca_dist, clusters = cluster_info)
    #umap_con <- clValid::connectivity(distance = umap_dist, clusters = cluster_info)
    #raw_con <- clValid::connectivity(distance = raw_dist, clusters = cluster_info)

    message("calculating Dunn index")
    #计算DUNN指数
    tsne_dunn <- clValid::dunn(distance = tsne_dist, clusters = as.numeric(as.factor(cluster_info)))
    pca_dunn <- clValid::dunn(distance = pca_dist, clusters = as.numeric(as.factor(cluster_info)))
    #umap_dunn <- clValid::dunn(distance = umap_dist, clusters = as.numeric(as.factor(cluster_info)))
    #raw_dunn <- clValid::dunn(distance = raw_dist, clusters = as.numeric(as.factor(cluster_info)))

    message("calculating Calinski Harabasz Index")
    #计算Calinski - Harabasz指数
    #raw_CHI <- fpc::calinhara(t(matrix),as.numeric(as.factor(cluster_info)))
    tsne_CHI <- fpc::calinhara(tsne_matrix,as.numeric(as.factor(cluster_info)))
    pca_CHI <- fpc::calinhara(pca_matrix,as.numeric(as.factor(cluster_info)))
    #umap_CHI <- fpc::calinhara(umap_matrix,as.numeric(as.factor(cluster_info)))

    message("calculating Davies-Bouldin Index")
    #计算戴维斯 - 布尔丁(DBI)指数
    #raw_DB <- clusterSim::index.DB(t(matrix), as.numeric(as.factor(cluster_info)))$DB
    tsne_DB <- clusterSim::index.DB(tsne_matrix, as.numeric(as.factor(cluster_info)))$DB
    pca_DB <- clusterSim::index.DB(pca_matrix, as.numeric(as.factor(cluster_info)))$DB
    #umap_DB <- clusterSim::index.DB(umap_matrix, as.numeric(as.factor(cluster_info)))$DB

    data_list <- list(
      ASW = list(
        PCA_asw = pca_average_silhouette,
        TSNE_asw = tsne_average_silhouette
      ),
      Dunn = list(
        PCA_Dunn = pca_dunn,
        TSNE_Dunn = tsne_dunn
      ),
      connectivity = list(
        PCA_connectivity = pca_con,
        TSNE_connectivity = tsne_con
      ),
      CHI = list(
        PCA_CHI = pca_CHI,
        TSNE_CHI = tsne_CHI
      ),
      DBI = list(
        PCA_DBI = pca_DB,
        TSNE_DBI = tsne_DB
      )
    )
    return(data_list)
}

#' Calculate the similarity in clustering evaluation metrics between the first dataset and other datasets
#'
#' @description
#' This method calculates the similarity between the clustering index of the first SCE object and the clustering indices of other SCE objects,
#' and normalizes the similarity to obtain the clustering similarity score between the first SCE object and other SCE objects.
#'
#' The calculated score represents the similarity in clustering performance between the original dataset and the simulated dataset,
#' with higher scores indicating greater similarity between the original and simulated data.
#'
#' @param SCEList SCEList is a list containing multiple SingleCellExperiment objects.Please ensure that each SCE object has a 'group' column in its colData.
#' @param TSNE_ntop A numerical scalar that specifies the number of features with the highest variance used for TSNE dimensionality reduction.
#' @param PCA_ntop A numerical scalar that specifies the number of features with the highest variance used for PCA dimensionality reduction.
#'
#' @return The value and list of clustering evaluation
#'
#' @import SingleCellExperiment
#'
#' @examples
#' data("ESCO_sce_293T_jurkat", "muscat_sce_293T_jurkat", "scDesign3_sce_293T_jurkat", "origin_sce_293T_jurkat")
#' SCEList <- list(origin=origin_sce_293T_jurkat, ESCO=ESCO_sce_293T_jurkat, muscat=muscat_sce_293T_jurkat, scDesign3 = scDesign3_sce_293T_jurkat)
#' test <- calculateClusterParameters(SCEList)
#'
#' @export
calculateClusterParameters <- function(SCEList, TSNE_ntop = 500, PCA_ntop = 500){

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

  if(nDatasets==1){
    stop("SCEList must be a named list, ",
         "must have at least two elements.", call. = FALSE)
  }


  result_list <- lapply(SCEList, function(sce){
    expr_data <- assays(sce)$counts
    group <- sce$group
    calculateCluster(expr_data, group, TSNE_ntop = TSNE_ntop, PCA_ntop = PCA_ntop)
  })

  score_list <- vector("list", length(result_list)-1)
  names(score_list) <- names(result_list)[-1]  # 保留原名称

  origin_params <- result_list[[1]]

  for (i in seq_along(score_list)) {
    # 提取当前元素的值
    current_pca_asw <- result_list[[i+1]]$ASW$PCA_asw
    current_tsne_asw <- result_list[[i+1]]$ASW$TSNE_asw
    #current_umap_asw <- result_list[[i+1]]$ASW$UMAP_asw

    current_pca_dunn <- result_list[[i+1]]$Dunn$PCA_Dunn
    current_tsne_dunn <- result_list[[i+1]]$Dunn$TSNE_Dunn

    current_pca_connectivity <- result_list[[i+1]]$connectivity$PCA_connectivity
    current_tsne_connectivity <- result_list[[i+1]]$connectivity$TSNE_connectivity

    current_pca_CHI <- result_list[[i+1]]$CHI$PCA_CHI
    current_tsne_CHI <- result_list[[i+1]]$CHI$TSNE_CHI

    current_pca_DBI <- result_list[[i+1]]$DBI$PCA_DBI
    current_tsne_DBI <- result_list[[i+1]]$DBI$TSNE_DBI

    # 计算差值
    diff_pca_asw <- current_pca_asw - origin_params$ASW$PCA_asw
    diff_tsne_asw <- current_tsne_asw - origin_params$ASW$TSNE_asw
    #diff_umap_asw <- current_umap_asw - origin_params$ASW$UMAP_asw

    diff_pca_dunn <- current_pca_dunn - origin_params$Dunn$PCA_Dunn
    diff_tsne_dunn <- current_tsne_dunn - origin_params$Dunn$TSNE_Dunn

    diff_pca_connectivity <- current_pca_connectivity - origin_params$connectivity$PCA_connectivity
    diff_tsne_connectivity <- current_tsne_connectivity - origin_params$connectivity$TSNE_connectivity

    diff_pca_chi <- current_pca_CHI - origin_params$CHI$PCA_CHI
    diff_tsne_chi <- current_tsne_CHI - origin_params$CHI$TSNE_CHI

    diff_pca_dbi <- current_pca_DBI - origin_params$DBI$PCA_DBI
    diff_tsne_dbi <- current_tsne_DBI - origin_params$DBI$TSNE_DBI

    # 保持原有结构存储结果
    score_list[[i]] <- list(
      ASW = list(
        PCA_asw = diff_pca_asw,
        TSNE_asw = diff_tsne_asw
        #UMAP_asw = diff_umap_asw
      ),
      Dunn = list(
        PCA_Dunn = diff_pca_dunn,
        TSNE_Dunn = diff_tsne_dunn
      ),
      connectivity = list(
        PCA_connectivity = diff_pca_connectivity,
        TSNE_connectivity = diff_tsne_connectivity
      ),
      CHI = list(
        PCA_CHI = diff_pca_chi,
        TSNE_CHI = diff_tsne_chi
      ),
      DBI = list(
        PCA_DBI = diff_pca_dbi,
        TSNE_DBI = diff_tsne_dbi
      )
    )
  }


  # 标准化与归一化函数
  scale_and_normalize <- function(data_matrix, center = TRUE, scale = TRUE) {

    if(all(data_matrix == data_matrix[1])){
      normalized_data=array(0, dim = length(data_matrix))
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
  pca_asw_matrix <- sapply(score_list, function(x) abs(x$ASW$PCA_asw))
  pca_asw_normalize <- scale_and_normalize(pca_asw_matrix)

  tsne_asw_matrix <- sapply(score_list, function(x) abs(x$ASW$TSNE_asw))
  tsne_asw_normalize <- scale_and_normalize(tsne_asw_matrix)

  pca_dunn_matrix <- sapply(score_list, function(x) abs(x$Dunn$PCA_Dunn))
  pca_dunn_normalize <- scale_and_normalize(pca_dunn_matrix)

  tsne_dunn_matrix <- sapply(score_list, function(x) abs(x$Dunn$TSNE_Dunn))
  tsne_dunn_normalize <- scale_and_normalize(tsne_dunn_matrix)

  pca_connectivity_matrix <- sapply(score_list, function(x) abs(x$connectivity$PCA_connectivity))
  pca_connectivity_normalize <- scale_and_normalize(pca_connectivity_matrix)

  tsne_connectivity_matrix <- sapply(score_list, function(x) abs(x$connectivity$TSNE_connectivity))
  tsne_connectivity_normalize <- scale_and_normalize(tsne_connectivity_matrix)

  pca_chi_matrix <- sapply(score_list, function(x) abs(x$CHI$PCA_CHI))
  pca_chi_normalize <- scale_and_normalize(pca_chi_matrix)

  tsne_chi_matrix <- sapply(score_list, function(x) abs(x$CHI$TSNE_CHI))
  tsne_chi_normalize <- scale_and_normalize(tsne_chi_matrix)

  pca_dbi_matrix <- sapply(score_list, function(x) abs(x$DBI$PCA_DBI))
  pca_dbi_normalize <- scale_and_normalize(pca_dbi_matrix)

  tsne_dbi_matrix <- sapply(score_list, function(x) abs(x$DBI$TSNE_DBI))
  tsne_dbi_normalize <- scale_and_normalize(tsne_dbi_matrix)


  normalized_list <- score_list

  #normalized_list
  for (i in seq_along(score_list)) {

    normalized_list[[i]]$ASW$PCA_asw <- pca_asw_normalize[[i]]
    normalized_list[[i]]$ASW$TSNE_asw <- tsne_asw_normalize[[i]]

    normalized_list[[i]]$Dunn$PCA_Dunn <- pca_dunn_normalize[[i]]
    normalized_list[[i]]$Dunn$TSNE_Dunn <- tsne_dunn_normalize[[i]]

    normalized_list[[i]]$connectivity$PCA_connectivity <- pca_connectivity_normalize[[i]]
    normalized_list[[i]]$connectivity$TSNE_connectivity <- tsne_connectivity_normalize[[i]]

    normalized_list[[i]]$CHI$PCA_CHI <- pca_chi_normalize[[i]]
    normalized_list[[i]]$CHI$TSNE_CHI <- tsne_chi_normalize[[i]]

    normalized_list[[i]]$DBI$PCA_DBI <- pca_dbi_normalize[[i]]
    normalized_list[[i]]$DBI$TSNE_DBI <- tsne_dbi_normalize[[i]]
  }

  for (i in seq_along(normalized_list)) {
    normalized_list[[i]]$score$PCA <- (1-normalized_list[[i]]$ASW$PCA_asw)+(1-normalized_list[[i]]$Dunn$PCA_Dunn)+
      (1-normalized_list[[i]]$connectivity$PCA_connectivity)+(1-normalized_list[[i]]$CHI$PCA_CHI)+
      (1-normalized_list[[i]]$DBI$PCA_DBI)

    normalized_list[[i]]$score$TSNE <- (1-normalized_list[[i]]$ASW$TSNE_asw)+(1-normalized_list[[i]]$Dunn$TSNE_Dunn)+
      (1-normalized_list[[i]]$connectivity$TSNE_connectivity)+(1-normalized_list[[i]]$CHI$TSNE_CHI)+
      (1-normalized_list[[i]]$DBI$TSNE_DBI)

    normalized_list[[i]]$score$sum <- normalized_list[[i]]$score$PCA + normalized_list[[i]]$score$TSNE

  }


  sum <- sapply(normalized_list, function(x) x$score$sum)

  # 构建 data.frame
  result_df <- data.frame(
    sum_score = sum,
    row.names = names(normalized_list)  # 设置行名为列表名
  )
  normalized_list$score <- result_df
  fin_list <- list(result_list = result_list, normalized_list = normalized_list)
  return(fin_list)


}


