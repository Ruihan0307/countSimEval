#' Construct data frame with pairwise statistics
#'
#' Construct a data frame containing statistics and p-values for pairwise
#' comparison of data sets.
#'
#' @param df The input data frame. Must contain at least a column named
#'   'dataset' and an additional column with values
#' @param column The name of the column(s) of \code{df} to be used as the basis
#'   for the comparison
#' @param permutationPvalues Whether or not to calculate p-values of statistics
#'   via permutation
#' @param nPermutations The number of permutations (only used if
#'   permutationPvalues = TRUE)
#' @param subsampleSize The number of observations for which certain
#'   (time-consuming) statistics will be calculated. The observations will be
#'   selected randomly among the rows of \code{df}
#' @param kmin,kfrac For statistics that require the extraction of k nearest
#'   neighbors of a given point, the number of neighbors will be max(kmin, kfrac
#'   * nrow(df))
#'
#' @return A data table with statistics and p-values for pairwise comparisons of
#'   data sets, based on the provided \code{column}
#'
#' @author Charlotte Soneson
#' @keywords internal
#'
#' @import dplyr
#' @importFrom utils combn
#'
makeDF <- function(df, column, permutationPvalues = FALSE, nPermutations = NULL,
                   subsampleSize, kmin, kfrac) {
  ## Generate data frame with comparisons of all pairs of data sets
  if (length(unique(df$dataset)) > 1) {
    ## Initialize data frame by populating it with all data set pairs
    ## ds_df <- data.frame(t(utils::combn(unique(df$dataset), 2)),
    ##                     stringsAsFactors = FALSE) %>%
    ##   dplyr::rename(dataset1 = X1, dataset2 = X2)


    unique_datasets <- unique(df$dataset)

    # 生成第一个数据集与其余所有数据集的成对组合
    ds_df <- expand.grid(dataset1 = unique_datasets[1],
                         dataset2 = unique_datasets[-1],
                         stringsAsFactors = FALSE)




    ## Calculate statistics for each pair of data sets
    ## 将n个放入的数据集分对，成对计算
    ds_res <- as.data.frame(t(vapply(seq_len(nrow(ds_df)), function(i) {
      dftmp <- df %>% dplyr::filter(dataset %in%
                                      c(ds_df$dataset1[i], ds_df$dataset2[i]))

      obs_stats <- calculateStats(df = dftmp, ds1 = ds_df$dataset1[i],
                                  ds2 = ds_df$dataset2[i], column = column,
                                  subsampleSize = subsampleSize,
                                  permute = FALSE, kmin = kmin, kfrac = kfrac,
                                  xmax = max(df[, column], na.rm = TRUE),
                                  xmin = min(df[, column], na.rm = TRUE))

      if (permutationPvalues) {
        perm_stats <- t(vapply(seq_len(nPermutations), function(s) {
          calculateStats(df = dftmp, ds1 = ds_df$dataset1[i],
                         ds2 = ds_df$dataset2[i], column = column,
                         subsampleSize = subsampleSize,
                         permute = TRUE, kmin = kmin, kfrac = kfrac,
                         xmax = max(df[, column], na.rm = TRUE),
                         xmin = min(df[, column], na.rm = TRUE))
        }, defaultStats(length(column), withP = FALSE)))
        ## Permutation p-values
        obs_stats["NNmismatchP"] <-
          signif(mean(c(obs_stats["NNmismatch"], perm_stats[, "NNmismatch"]) >=
                        obs_stats["NNmismatch"]), 3)
        obs_stats["avesilhP"] <-
          signif(mean(c(obs_stats["avesilh"], perm_stats[, "avesilh"]) >=
                        obs_stats["avesilh"]), 3)
        obs_stats["avesilhlocalP"] <-
          signif(mean(c(obs_stats["avesilhlocal"],
                        perm_stats[, "avesilhlocal"]) >=
                        obs_stats["avesilhlocal"]), 3)
        if (length(column) == 1)
          obs_stats["diffareaP"] <-
          signif(mean(c(obs_stats["diffarea"], perm_stats[, "diffarea"]) >=
                        obs_stats["diffarea"]), 3)
      }
      obs_stats
    }, defaultStats(length(column), withP = permutationPvalues))))

    ## Populate the output table



    if (length(column) == 1) {
      ds_df$`K-S statistic` <- ds_res$ksstatistic
      ds_df$`K-S p-value` <- ds_res$kspvalue
      ds_df$`Scaled area between eCDFs` <- ds_res$diffarea
      if (permutationPvalues) {
        ds_df$`Scaled area between eCDFs perm p-value` <- ds_res$diffareaP
      }
      ds_df$`Runs statistic` <- ds_res$runsstatistic
      ds_df$`Runs p-value` <- ds_res$runspvalue
    }
    ds_df$`NN rejection fraction` <- ds_res$NNmismatch
    if (permutationPvalues) {
      ds_df$`NN rejection fraction perm p-value` <- ds_res$NNmismatchP
    }
    ds_df$`Average silhouette width` <- ds_res$avesilh
    if (permutationPvalues) {
      ds_df$`Average silhouette width perm p-value` <- ds_res$avesilhP
    }
    ds_df$`Average local silhouette width` <- ds_res$avesilhlocal
    if (permutationPvalues) {
      ds_df$`Average local silhouette width perm p-value` <-
        ds_res$avesilhlocalP
    }
    ds_df$`Multivariate KS statistic` <- ds_res$MultiKSstatistic
    if (permutationPvalues) {
      ds_df$`Multivariate KS p-value` <- ds_res$MultiKSpvalue
    }
    ds_df$`KDE T-statistic` <- ds_res$KDETstat
    ds_df$`KDE Z-statistic` <- ds_res$KDEZstat
    if (permutationPvalues) {
      ds_df$`KDE p-value` <- ds_res$KDEpvalue
    }
    ## Return output table
    ##DT::datatable(ds_df, options = list(
    ##  scrollX = TRUE,
    ##  columnDefs = list(list(className = "dt-center", targets = 3:ncol(ds_df)))
    ##))
    return(ds_df)
  }
}

#' Define table descriptions
#'
#' Generate the text that describes the content of the tables generated by
#' \code{makeDF}.
#'
#' @param calculateStatistics Whether or not statistics and p-values are
#'   calculated
#' @param subsampleSize The number of observations for which certain
#'   (time-consuming) statistics will be calculated
#' @param kmin,kfrac For statistics that require the extraction of k nearest
#'   neighbors of a given point, the number of neighbors will be max(kmin, kfrac
#'   * nrow(df))
#' @param obstype The type of observation (e.g., sample, feature, sample pair)
#' @param aspect The name of the aspect of interest
#' @param minvalue,maxvalue The minimal and maximal value of the aspect of
#'   interest, used for scaling of the x axis when calculating the area between
#'   the eCDFs
#' @param permutationPvalues Whether or not to calculate p-values of statistics
#'   via permutation
#' @param nPermutations The number of permutations (only used if
#'   permutationPvalues = TRUE)
#' @param nDatasets The number of data sets that are being compared
#'
#' @return A list with two text strings in markdown format: one for tables based
#'   on a single data column, and one for tables based on two data columns
#'
#' @keywords internal
#' @author Charlotte Soneson
#'
defineTableDesc <- function(calculateStatistics, subsampleSize, kfrac, kmin,
                            obstype, aspect, minvalue, maxvalue,
                            permutationPvalues, nPermutations, nDatasets) {
  if (nDatasets < 2) {
    tabledesc <- tabledesc2d <-
      "No statistics were calculated, since there is only one data set."
  } else {
    if (calculateStatistics) {
      if (permutationPvalues) {
        permtext <- paste0("The permutation p-values below were estimated ",
                           "based on ", nPermutations, " permutations.")
      } else {
        permtext <- ""
      }
      tabledesc <- paste0("The table below contains a range of quantitative ",
                          "statistics and test results evaluating the degree ",
                          "of similarity between each pair of data sets based ",
                          "on the ", aspect, ". The following statistics and ",
                          "test results are included:

- *K-S statistic/K-S p-value*: the [Kolmogorov-Smirnov statistic]",
"(https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test) ",
"[@Kolmogorov1933, @Smirnov1948], i.e., the maximal distance between the two ",
"empirical cumulative distribution functions (eCDFs), and the associated ",
"p-value.
- *Scaled area between eCDFs*: the absolute area between the eCDFs for the ",
"two data sets. This is calculated by first subtracting one eCDF from the ",
"other and taking the absolute value of the differences, followed by scaling ",
"the x-axis so that the difference between the largest and smallest value ",
"across all data sets is equal to 1 (here, subtracting ", minvalue, " and ",
"dividing by ", maxvalue - minvalue, "), and finally calculating the area ",
"under the resulting non-negative curve. This gives a scaled area with ",
"values in the range [0, 1], where large values indicate large differences ",
"between the ", aspect, " distributions.
- *Runs statistic/Runs p-value*: the statistic and p-value from a one-sided ",
"[Wald-Wolfowitz runs test]",
"(https://en.wikipedia.org/wiki/Wald%E2%80%93Wolfowitz_runs_test) ",
"[@WaldWolfowitz1940]. The values from both data sets are first sorted in ",
"increasing order and converted into a binary sequence based on whether or ",
"not they are derived from the first data set. The runs test is applied to ",
"this sequence. The test is one-sided, with the alternative hypothesis being ",
"that there is a trend towards values from the same data set occurring next ",
"to each other in the sequence.
- *NN rejection fraction*: this statistic is similar to one used by the ",
"[kBET](https://github.com/theislab/kBET) R package. First, a subset of R ",
obstype, "s are selected (where R is chosen to be the smallest of ",
subsampleSize, " and the total number of ", obstype, "s in the two data sets)",
", and for each of these ", obstype, "s the k nearest neighbors are found ",
"(where k is chosen to be ", kfrac * 100, "% of the total number of ", obstype,
"s from the two data sets, however at least ", kmin, "). Then, a chi-square ",
"test is used to evaluate whether the composition of data sets from which ",
"the k nearest neighbors of each selected ", obstype, "s comes differs ",
"significantly from the overall composition of ", obstype, "s from the two ",
"data sets. The NN rejection fraction is the fraction of selected ", obstype,
"s for which the chi-square test rejects the null hypothesis at a ",
"significance level of 5%.
- *Average silhouette width*: as for the NN rejection fraction, a subset of R ",
obstype, "s are selected. For each of these, the (Euclidean) distances to all ",
"other ", obstype, "s are calculated, and the [silhouette width]",
"(https://en.wikipedia.org/wiki/Silhouette_(clustering)) ",
"[@Rousseeuw1987silhouette] is calculated based on these distances. The final",
" statistic is the average silhouette width for the R ", obstype, "s.
- *Average local silhouette width*: similar to the average silhouette width, ",
"but for each ", obstype, ", only the k' closest ", obstype, "s from a data ",
"set are used to calculate the silhouette widths. Here, k' for data set i is ",
"defined as k * p_i, where k is as above and p_i is the fraction of the total ",
"number of ", obstype, "s that come from data set i.
- For the scaled area between eCDFs, NN rejection fraction and average ",
"silhouette widths, permutation p-values can be calculated by permuting the ",
"data set labels. These values are only available if the 'permutationPvalues' ",
"argument to 'countsimQCReport' was set to TRUE. ", permtext,"

It should be noted that the reported statistics in some cases are highly ",
"dependent on the number of underlying observations, and especially with a ",
"large number of observations, very small p-values can correspond to small ",
"effective differences. Thus, interpretation should ideally be guided by both ",
"quantitative and qualitative, visual information.")
      tabledesc2d <- paste0("The table below contains a range of quantitative ",
"statistics and test results evaluating the degree of similarity between each ",
"pair of data sets based on the ", aspect, ". The following statistics and ",
"test results are included:

- *NN rejection fraction*: this statistic is similar to one used by the ",
"[kBET](https://github.com/theislab/kBET) R package. First, a subset of R ",
obstype, "s are selected (where R is chosen to be the smallest of ",
subsampleSize, " and the total number of ", obstype, "s in the two data sets)",
", and for each of these ", obstype, "s the k nearest neighbors are found ",
"(where k is chosen to be ", kfrac * 100, "% of the total number of ", obstype,
"s from the two data sets, however at least ", kmin, "). Then, a chi-square ",
"test is used to evaluate whether the composition of data sets from which the ",
"k nearest neighbors of each selected ", obstype, "s comes differs ",
"significantly from the overall composition of ", obstype, "s from the two ",
"data sets. The NN rejection fraction is the fraction of selected ", obstype,
"s for which the chi-square test rejects the null hypothesis at a ",
"significance level of 5%.
- *Average silhouette width*: as for the NN rejection fraction, a subset of R ",
obstype, "s are selected. For each of these, the (Euclidean) distances to ",
"all other ", obstype, "s are calculated, and the [silhouette width]",
"(https://en.wikipedia.org/wiki/Silhouette_(clustering)) ",
"[@Rousseeuw1987silhouette] is calculated based on these distances. The final ",
"statistic is the average silhouette width for the R ", obstype, "s.
- *Average local silhouette width*: similar to the average silhouette width, ",
"but for each ", obstype, ", only the k' closest ", obstype, "s from a data ",
"set are used to calculate the silhouette widths. Here, k' for data set i is ",
"defined as k * p_i, where k is as above and p_i is the fraction of the ",
"total number of ", obstype, "s that come from data set i.
- For the NN rejection fraction and average silhouette widths, permutation ",
"p-values can be calculated by permuting the data set labels. These values ",
"are only available if the 'permutationPvalues' argument to ",
"'countsimQCReport' was set to TRUE. ", permtext)
    } else {
      tabledesc <- tabledesc2d <-
        paste0("No statistics were calculated, since the ",
               "'calculateStatistics' argument to 'countsimQCReport()' was ",
               "set to FALSE. To perform pairwise quantitative comparisons ",
               "between data sets, set this argument to TRUE. Note, however, ",
               "that this increases the runtime significantly.")
    }
  }

  list(tabledesc = tabledesc, tabledesc2d = tabledesc2d)
}
