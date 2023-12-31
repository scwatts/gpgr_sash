% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmd.R
\name{cancer_rmd}
\alias{cancer_rmd}
\title{Generate UMCCR Cancer Report}
\usage{
cancer_rmd(
  af_global,
  af_keygenes,
  batch_name,
  conda_list,
  img_dir,
  key_genes,
  oncoviral_breakpoints_tsv,
  oncoviral_present_viruses,
  purple_germ_cnv,
  purple_purity,
  purple_qc,
  purple_som_cnv,
  purple_som_gene_cnv,
  purple_som_snv_vcf,
  somatic_snv_summary,
  somatic_snv_vcf,
  somatic_sv_tsv,
  somatic_sv_vcf,
  result_outdir,
  tumor_name,
  out_file = NULL,
  quiet = FALSE
)
}
\arguments{
\item{af_global}{Path to \code{af_tumor.txt} file.}

\item{af_keygenes}{Path to \code{af_tumor_keygenes.txt} file.}

\item{batch_name}{Name of batch sample.}

\item{conda_list}{Path to \code{conda_pkg_list.txt} file.}

\item{img_dir}{Path to directory containing PURPLE plots.}

\item{key_genes}{Path to UMCCR cancer gene file.}

\item{oncoviral_breakpoints_tsv}{Path to \code{oncoviruses/oncoviral_breakpoints.tsv}.}

\item{oncoviral_present_viruses}{Path to \code{oncoviruses/present_viruses.txt}.}

\item{purple_germ_cnv}{Path to \code{purple.cnv.germline.tsv}.}

\item{purple_purity}{Path to \code{purple.purity.tsv}.}

\item{purple_qc}{Path to \code{purple.qc}.}

\item{purple_som_cnv}{Path to \code{purple.cnv.somatic.tsv}.}

\item{purple_som_gene_cnv}{Path to \code{purple.cnv.gene.tsv}.}

\item{purple_som_snv_vcf}{Path to \code{purple.somatic.vcf.gz}.}

\item{somatic_snv_summary}{Path to \code{somatic_snv_summary.json} JSON.}

\item{somatic_snv_vcf}{Path to \code{somatic-PASS.vcf.gz} SNV VCF.}

\item{somatic_sv_tsv}{Path to \code{manta.tsv} TSV file.}

\item{somatic_sv_vcf}{Path to \code{manta.vcf.gz} VCF file.}

\item{result_outdir}{Path to directory to write tidy JSON/TSV results.}

\item{tumor_name}{Name of tumor sample.}

\item{out_file}{Path to output HTML file (needs '.html' suffix) (def: \verb{\{tumor_name\}_cancer_report.html}).}

\item{quiet}{Suppress log printing during rendering.}
}
\value{
Path to rendered HTML report.
}
\description{
Generates a UMCCR Cancer Report. It does so with the following steps:
\enumerate{
\item move the img_dir into 'tmp/img_dir'
\item copy the rmd into 'tmp/cancer_report.Rmd'
\item render the rmd inside 'tmp/'
\item return the path to the output HTML
}
}
