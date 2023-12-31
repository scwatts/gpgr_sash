% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sv.R
\name{process_sv}
\alias{process_sv}
\title{Process Structural Variants}
\usage{
process_sv(x)
}
\arguments{
\item{x}{Path to \code{manta.tsv} output by umccrise.}
}
\value{
A list with melted/unmelted tibbles (these are NULL if TSV file was empty).
}
\description{
Processes the Manta TSV file output by umccrise.
}
\examples{
x <- system.file("extdata/umccrise/sv/manta.tsv", package = "gpgr")
(sv <- process_sv(x))
}
