% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotLog2fcHeatmap.R
\name{plotLog2fcHeatmap}
\alias{plotLog2fcHeatmap}
\title{Generate a heatmap and average profile plot of log2 fold changes around summit regions}
\usage{
plotLog2fcHeatmap(
  tx_bam,
  in_bam,
  summit_file_paths,
  peak_file_paths,
  plotting_range = 4001,
  window_size = 150,
  step_size = 5,
  core_count = 8,
  chromosomes = paste0("chr", 1:22),
  color_breaks = c(-0.5, 0, 0.5, 1, 2)
)
}
\arguments{
\item{tx_bam}{Character. Path to treatment BAM file.}

\item{in_bam}{Character. Path to input/control BAM file.}

\item{summit_file_paths}{Named character vector. Paths to summit BED files for each group (e.g., TRESLIN, MTBP, Both).}

\item{peak_file_paths}{Named character vector. Paths to reproducible peak BED files for each group.}

\item{plotting_range}{Integer. Width (in bp) around each summit for signal profiling. Default: 4001.}

\item{window_size}{Integer. Size of the rolling window in bp. Default: 150.}

\item{step_size}{Integer. Step size for sliding window. Default: 5.}

\item{core_count}{Integer. Number of CPU cores for parallel processing. Default: 8.}

\item{chromosomes}{Character vector. Chromosomes to include. Default: paste0("chr", 1:22).}

\item{color_breaks}{Numeric vector. Values to define color scale (length should match number of colors). Default: c(-0.5, 0, 0.5, 1, 2).}
}
\value{
A list with two elements:
\describe{
\item{average_profile_plot}{A ggplot2 line plot showing average log2 fold change per position per group.}
\item{heatmap_plot}{A ComplexHeatmap object showing per-region log2 fold change.}
}
}
\description{
This function computes and visualizes signal profiles (log2 fold change) around summit regions.
It outputs a ComplexHeatmap object and a ggplot2 line plot showing average signal per group.
}
