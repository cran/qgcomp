#' @title Well water data
#' 
#' @description Simulated well water measurements in North Carolina: 16 metals, 6 water chemistry 
#' measures, and 2 health outcomes (y = continuous; disease_state = binary/time-to-event
#' in combination with disease_time)
#'
#' A dataset containing well water measurements and health outcomes for 253 individuals. 
#' All continuous variables are standardized to have mean 0, standard deviation 1.
#' @usage data(metals)
#' @format A data frame with 253 rows and 24 variables:
#' \describe{
#'   \item{y}{continuous birth outcome}
#'   \item{disease_state}{binary outcome}
#'   \item{disease_time}{time-to-disease_state: survival outcome censored at approximately the median}
#'   \item{arsenic}{metal}
#'   \item{barium}{metal}
#'   \item{cadmium}{metal}
#'   \item{calcium}{metal}
#'   \item{chloride}{metal}
#'   \item{chromium}{metal}
#'   \item{copper}{metal}
#'   \item{iron}{metal}
#'   \item{lead}{metal}
#'   \item{magnesium}{metal}
#'   \item{manganese}{metal}
#'   \item{mercury}{metal}
#'   \item{selenium}{metal}
#'   \item{silver}{metal}
#'   \item{sodium}{metal}
#'   \item{zinc}{metal}
#'   \item{mage35}{Binary covariate: maternal age > 35}
#'   \item{nitrate}{water chemistry measure}
#'   \item{nitrite}{water chemistry measure}
#'   \item{sulfate}{water chemistry measure}
#'   \item{ph}{water chemistry measure}
#'   \item{total_alkalinity}{water chemistry measure}
#'   \item{total_hardness}{water chemistry measure}
#' }
"metals"
