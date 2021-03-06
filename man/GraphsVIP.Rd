% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GraphsVIP.R
\name{GraphsVIP}
\alias{GraphsVIP}
\title{Variable importance plots (VIP)}
\usage{
GraphsVIP(data, name, groupnames, tsf = "clr", top = 30, QCs = FALSE)
}
\arguments{
\item{data}{Data table with variables (metabolites) in columns. Samples in rows are sorted according to specific groups.}

\item{name}{A character string or expression indicating a name of data set. It occurs in names of every output.}

\item{groupnames}{A character vector defining specific groups in data. Every string must be specific for each group and they must not overlap.}

\item{tsf}{Data transformation must be defined by "clr" (default), "log", "log10", "PQN", "lnPQN", "pareto" or "none". See Details.}

\item{top}{How many most important variables should be in zoomed VIP plot and biplot? The default is 30.}

\item{QCs}{logical. If FALSE (default) quality control samples (QCs) are automatically distinguished and skipped.}
}
\value{
VIP plots, score plots and biplots of PLS-DA.

Excel file with mean, sd and angles of rays from group center for every variable.
}
\description{
Makes variable importance plots (VIP) of partial least squares - discriminant analysis (PLS-DA), also displays score plots and biplots.
}
\details{
Data transformation: with "clr" clr trasformation is used (see References), with "log" natural logarithm is used, with "log10" decadic logarithm is used, with "pareto" data are only scaled by Pareto scaling, with "PQN" probabilistic quotient normalization is done, with "lnPQN" natural logarithm of PQN transformed data is done, with "none" no tranformation is done.

VIP plot can be used only for comparison of two groups. If there is more groups in data, all possible combinations of pairs are evaluated.

Up to twenty different groups can be distinguished in data (including QCs).
}
\examples{
data=metabol
name="Metabolomics"    #name of the project
groupnames=c("Con","Pat","QC")
GraphsVIP(data,name,groupnames)
}
\references{
Aitchison, J. (1986) The Statistical Analysis of Compositional Data Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). p. 416.

Wold, S. et al. (2001) PLS-regression: a basic tool of chemometrics, CHEMOLAB 58 (2), p. 109-130.

Chong, I.G., Jun, C.H. (2005) Performance of some variable selection methods when multicollinearity is present, CHEMOLAB 78, p. 103-112.
}
