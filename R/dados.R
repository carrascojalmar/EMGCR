#' Liver Cancer Data
#'
#' A sample of 2,766 patients diagnosed with liver cancer between 2012 and 2016, whose cancer grades were well identified. Available individual-level covariates include age at diagnosis, sex, histological grade of the liver cancer, number of relapses, and median household income.
#'
#' The grade of the disease is categorized into four levels:
#' \enumerate{
#'   \item Grade I – Well differentiated
#'   \item Grade II – Moderately differentiated
#'   \item Grade III – Poorly differentiated
#'   \item Grade IV – Undifferentiated/anaplastic
#' }
#'
#' @docType data
#' @name liver
#' @usage data(liver)
#' @format A data frame with 2,766 observations and 10 variables:
#' \describe{
#'   \item{ID}{Unique patient identifier}
#'   \item{age}{Age as a factor (e.g., age group)}
#'   \item{ageNumeric}{Age as a numeric variable}
#'   \item{grade}{Histologic grade of liver cancer (I, II, III, IV)}
#'   \item{medh}{Median household income as factor (possibly grouped)}
#'   \item{medhNumeric}{Median household income as numeric}
#'   \item{relapse}{Number of relapses after first diagnosis}
#'   \item{sex}{Sex of patient: 1 = male, 0 = female}
#'   \item{status}{Event indicator: 1 = death, 0 = censored}
#'   \item{time}{Survival time in months}
#' }
#'
#' @examples
#' data(liver)
#' head(liver)
NULL
