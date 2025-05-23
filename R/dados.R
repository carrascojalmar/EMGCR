#' Liver cancer data
#'
#' @description A sample of 2,766 patients who have been diagnosed with liver cancer between 2012 and 2016 whose cancer grades are well identified. Available individual-level covariates include age at diagnosis, sex, grade of liver cancer, count of relapse, and median household income. The grade of the disease is categorized into four levels: well-differentiated (Grade I), moderately differentiated (Grade II), poorly differentiated (Grade III), and undifferentiated/anaplastic (Grade IV).
#'
#' @docType data
#' @name liver
#' @usage data(liver)
#' @format This data frame contains the following columns:
#' \itemize{
#' \item {time:} {survival time in month}
#' \item {status:} {censored = 0, dead = 1}
#' \item {sex:} {1 if male, 0 if female}
#' \item {age:} {Age at diagnosis }
#' \item {meth:} {Median household income of the subject }
#' \item {relapse:} {Count of relapses after first diagnosis }
#' \item {grade:} {Histologic grade of liver cancer}
#' }

#' @examples
#' data(liver)
#' head(liver)
#'
NULL
