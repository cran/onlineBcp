#' GC-corrected data for copy number variation
#'
#' A dataset containing the raw data and GC-corrected/normalized data
#'
#' @format A data frame with 14189 rows and 2 variables:
#' \describe{
#'   \item{raw.count}{raw read counts}
#'   \item{normalized.count}{normalized read counts}
#' }
"cnv_H2347"
#' Transformed aCGH data
#'
#' A dataset containing the tranformed aCGH data from the genome of the fibroblast cell line GM02948
#'
#' @format A data frame with 2046 rows and 1 variable:
#' \describe{
#' \item{transNorm}{normalized aCGH intensity}
#' }
"aCGH"
#' US COVID-19 data
#'
#' A dataset containing new daily cases in the United States downloaded from the
#' World Health Organization on August 25, 2020
#'
#' @format A data frame with 219 rows and 8 variables
#' \describe{
#' \item{Date_reported}{The report date}
#' \item{Country_code}{The code for country}
#' \item{Country}{Country in full name}
#' \item{WHO_region}{Geographic region defined by WHO}
#' \item{New_cases}{New COVID-19 cases}
#' \item{Cumulative_cases}{Cumulative COVID-19 cases}
#' \item{New_deaths}{New COVID-19 deaths}
#' \item{Cumulative_deaths}{Cumulative COVID-19 deaths}
#' }
"covid"
