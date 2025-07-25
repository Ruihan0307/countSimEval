% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcDataProperty.R
\name{calcDataProperty}
\alias{calcDataProperty}
\title{Calculate data attribute statistics}
\usage{
calcDataProperty(matList, maxNForDisp, maxNForCorr, subsampleSize, kfrac, kmin)
}
\arguments{
\item{matList}{List of gene expression count matrices}

\item{maxNForDisp}{The maximal number of samples that will be used to estimate dispersions. By default, all samples are used. This can be lowered to speed up calculations (and obtain approximate results) for large data sets.}

\item{maxNForCorr}{The maximal number of samples (features) for which pairwise correlation coefficients will be calculated. If the number of samples (features) exceeds this number, they will be randomly subsampled.}

\item{subsampleSize}{The number of randomly selected observations (samples, features or pairs of samples or features) for which certain (time-consuming) statistics will be calculated.}

\item{kfrac}{For statistics that require the extraction of the k nearest neighbors of a given point, the number of neighbors will be max(kmin, kfrac * nrow(df))}

\item{kmin}{For statistics that require the extraction of the k nearest neighbors of a given point, the number of neighbors will be max(kmin, kfrac * nrow(df))}
}
\value{
Calculate data attribute statistics
}
\description{
Calculates the similarity of data attributes between the first matrix and the other matrix.
}
