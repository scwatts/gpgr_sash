% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/purple.R
\name{purple_kataegis}
\alias{purple_kataegis}
\title{Get PURPLE Kataegis Regions}
\usage{
purple_kataegis(x)
}
\arguments{
\item{x}{Path to the \code{purple.somatic.vcf.gz} file.}
}
\value{
A list with a tibble containing variants in kataegis clusters and
various metrics for each variant, and a tibble with a description of each
metric.
}
\description{
Reads the \code{purple.somatic.vcf.gz} file and extracts variants
within kataegis regions.
}
\examples{
x <- system.file("extdata/purple/purple.somatic.vcf.gz", package = "gpgr")
(k <- purple_kataegis(x))
}
