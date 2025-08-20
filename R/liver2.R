#' Liver Cancer Data 2
#'
#' A sample of 1,736 patients diagnosed with liver cancer between 2012 and 2016, whose cancer grades were well identified. Available individual-level covariates include age at diagnosis, sex, histological grade of the liver cancer, number of relapses, and median household income.
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
#' @name liver2
#' @usage data(liver2)
#' @format A data frame with 1,736 observations and 10 variables:
#' \describe{
#'   \item{ID}{Unique patient identifier}
#'   \item{time}{Survival time in months}
#'   \item{status}{Event indicator: 1 = death, 0 = censored}
#'   \item{sex}{Sex of patient: 1 = male, 0 = female}
#'   \item{age}{Age as a factor (e.g., age group)}
#'   \item{medh}{Median household income as factor (possibly grouped)}
#'   \item{grade}{Histologic grade of liver cancer (I, II, III, IV)}
#'   \item{chemo}{chemotherap: 1 = yes}
#'   \item{radio}{radiation: 1 = yes}



#' }
#'
#' @examples
#' data(liver2)
#' head(liver2)
NULL
