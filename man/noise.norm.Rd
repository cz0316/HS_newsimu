% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/1004new_noise_simu.R
\name{noise.norm}
\alias{noise.norm}
\title{Title swap noise function 2}
\usage{
noise.norm(
  noise = 0.2,
  cut = 5,
  lambda = 0.7,
  spots,
  se,
  ns,
  se.p,
  ns.p,
  se.sd = 0.2,
  se.mu,
  ns.sd = 0.2,
  ns.mu,
  ptn,
  png_dir
)
}
\arguments{
\item{noise}{swap porpotion.}

\item{cut}{cur ,the threshold.}

\item{lambda}{The probability of the sample() function process.}

\item{spots}{The number of spots.}

\item{se}{The number of  spatially variable genes (SVGs).}

\item{ns}{The number of  non-SVGs.}

\item{se.p}{A number. The zero proportion of zero generation process of
SVGs in the streak area.}

\item{ns.p}{A number. The zero proportion of zero generation process of
non-SVGs and SVGs in the non-streak area.}

\item{se.mu}{For SVGs, the lambda para in the poisson distribution
or the mu para in the NB distribution.}

\item{ns.mu}{For non-SVGs and SVGs in the non-streak area, the lambda para in the poisson distribution
or the mu para in the NB distribution.}

\item{ptn}{The file name of the pattern png.}

\item{png_dir}{The dir of the png files.}

\item{type}{A character.'ZINB' or 'ZIP. The default value is ZINB.}

\item{se.size}{The size para of SVGs in the NB distribution.}

\item{ns.size}{For non-SVGs and SVGs in the non-streak area,
the size para in the NB distribution.}
}
\value{
A data.frame of gene expression counts and coordinates.
}
\description{
The swap noise function2, swap noise, in the  logNormal with cut.
}
