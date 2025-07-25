#' Calculate the number of DEGs between groups
#'
#' @param expr_data Gene expression matrix
#' @param group cell group
#' @param P.Value P.Value threshold for screening differentially expressed genes
#' @param logFC LogFC threshold for screening differentially expressed genes
#'
#' @return the number of DEGs
#'
#'
#' @export
calculateDEGs <- function(expr_data, group, P.Value, logFC){

  message(paste0("There are ",nlevels(factor(group))," groups"))

  dgelist <- DGEList(counts = expr_data, group = group)
  # 归一化，TMM 方法
  dgelist <- calcNormFactors(dgelist,method = "TMM")

  # 注意：归一化并不会直接在counts数值上修改，而是归一化系数会被自动存在 d$samples$norm.factors

  # 顺便介绍一下归一化的意义
  # 归一化不是绝对必要的，但是推荐进行归一化。
  # 有重复的样本中，应该不具备生物学意义的外部因素会影响单个样品的表达
  # 例如中第一批制备的样品会总体上表达高于第二批制备的样品，
  # 假设所有样品表达值的范围和分布都应当相似，
  # 需要进行归一化来确保整个实验中每个样本的表达分布都相似。


  # 创建设计矩阵，用于指定差异分析模型
  design <- model.matrix(~0 + factor(group))
  rownames(design) <- colnames(dgelist)
  colnames(design) <- levels(factor(group))

  # 估计数据的离散度 —— common离散度、trended离散度、tagwise离散度
  dgelist <- estimateDisp(dgelist, design, robust=T)
  #plotBCV(dgelist)
  ## To perform likelihood ratio test：scRNA-seq and no replicates data
  fit <- glmFit(dgelist, design, robust=T)

  n <- nlevels(factor(group))
  contrast <- list()
  for (i in 1:(n - 1)) {
    # 内层循环，遍历从 i+1 到 n 的组
    for (j in (i + 1):n) {
      # 生成一个长度为 n 的全 0 向量
      vec <- rep(0, n)
      # 第 i 个位置设为 -1
      vec[i] <- -1
      # 第 j 个位置设为 1
      vec[j] <- 1
      # 生成对比名称
      contrast_name <- paste0("Gp", i, "vsGp", j)
      # 将对比向量添加到列表中
      contrast[[contrast_name]] <- vec
    }
  }

  data_list <- lapply(contrast, function(contrast) {
    lt <- glmLRT(fit, contrast=contrast)
    tempDEG <- topTags(lt, n = nrow(dgelist$counts))

    origin_DEGs <- as.data.frame(tempDEG)

    #对up和down基因进行选择
    k1 <- (origin_DEGs$PValue < P.Value) & (origin_DEGs$logFC < -logFC)
    k2 <- (origin_DEGs$PValue < P.Value) & (origin_DEGs$logFC > logFC)
    origin_DEGs <- mutate(origin_DEGs, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
    return(origin_DEGs)

  })
}

#' Calculate the difference in DEGs between the first dataset and other datasets
#'
#' @description
#' This method calculates the differentially expressed genes of different groups between the first SCE object and each other SCE object.
#' The calculated score represents the difference in the number of differentially expressed genes between the original dataset and the simulated dataset.
#' The lower the score, the more similar the original and simulated data are
#'
#' @param P.Value P.Value threshold for screening differentially expressed genes
#' @param logFC LogFC threshold for screening differentially expressed genes
#' @param SCEList SCEList is a list containing multiple SingleCellExperiment objects. Please ensure that each SCE object has a 'group' column in its colData, and that the number of types in each SCE object group is the same.
#'
#' @return List of DEG quantities
#' @import SingleCellExperiment
#'
#' @examples
#' data("ESCO_sce_293T_jurkat", "muscat_sce_293T_jurkat", "scDesign3_sce_293T_jurkat", "origin_sce_293T_jurkat")
#' SCEList <- list(origin=origin_sce_293T_jurkat, ESCO=ESCO_sce_293T_jurkat, muscat=muscat_sce_293T_jurkat, scDesign3 = scDesign3_sce_293T_jurkat)
#' test <- calculateDEGsParameters(SCEList)
#'
#' @export

calculateDEGsParameters <- function(SCEList, P.Value = 0.05, logFC = 1){

  ##数据集数
  nDatasets <- length(SCEList)
  message(paste("There are a total of ", nDatasets, "datasets"))

  ## sce对象必须有counts和group
  DEGS <- lapply(SCEList, function(sce){
    expr_data <- assays(sce)$counts
    group <- sce$group
    calculateDEGs(expr_data, group, P.Value = P.Value, logFC = logFC)
  })
  new_list <- setNames(vector("list", length(DEGS) - 1), names(DEGS)[-1])

  for(name in names(DEGS[[1]])){
    up_count <- sum(DEGS[[1]][[name]]$change == "up")
    down_count <- sum(DEGS[[1]][[name]]$change == "down")

    for (i in seq_along(new_list)) {
      new_list[[i]][[name]]$up = abs((sum(DEGS[[i+1]][[name]]$change == "up")-up_count)*2)/(sum(DEGS[[i+1]][[name]]$change == "up")+up_count)
      new_list[[i]][[name]]$down = abs((sum(DEGS[[i+1]][[name]]$change == "down")-down_count)*2)/(sum(DEGS[[i+1]][[name]]$change == "down")+down_count)
    }
  }
  new_list <- lapply(new_list, function(name){
    total_up <- 0
    total_down <- 0
    for(vs in name){
      total_up <- vs$up+total_up
      total_down <- vs$down+total_down
    }
    name$score$total_up <- total_up
    name$score$total_down <- total_down
    return(name)
  })

  total_up <- sapply(new_list, function(x) x$score$total_up)
  total_down <- sapply(new_list, function(x) x$score$total_down)

  # 构建 data.frame
  result_df <- data.frame(
    total_up = total_up,
    total_down = total_down,
    row.names = names(new_list)  # 设置行名为列表名
  )
  new_list$score <- result_df
  result_list <- list(DEGs = DEGS, results = new_list)
  return(result_list)
}
