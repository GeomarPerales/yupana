#' Log gumbel distribution function for 24-hour maximum precipitation
#'
#' Log gumbel distribution function for 24-hour maximum precipitation for the study
#' of extreme precipitation events.
#'
#' @param x Vector or dataframe of daily maximum precipitation values
#' @param tr OPTIONAL, Vector of return period values, default parameter is NULL.
#'  Typical values for TR are 2, 5, 10, 20, 50, 100, 200, 500, 1000 years.
#'
#' @export
#' @importFrom stats na.omit sd
#'
#' @author Geomar Perales Apaico
#'
#' @name dislgum

dislgum <- function(x, tr = NULL) {
  if (is.data.frame(x)) {
    if (ncol(x) == 2) {
      x <- as.numeric(x[, 2])
    } else if (ncol(x) == 1) {
      x <- as.numeric(x)
    } else {
      stop("values not defined")
    }
  } else if (is.numeric(x)) {
    x <- x
  } else {
    stop("values not defined")
  }

  x <- na.omit(x)
  xlog <- log(x)
  mu_xlog <- mean(xlog)
  sigma_xlog <- sd(xlog)

  gamma_euler <- 0.5772156649

  alfa <- sigma_xlog * sqrt(6) / pi
  beta <- mu_xlog - gamma_euler * alfa

  TR <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)

  P <- 1 - (1 / TR)

  Y_gt <- -log(-log(P))
  Q_gumbel_xlog <- beta + alfa * Y_gt
  Q_loggumbel <- exp(Q_gumbel_xlog)

  if (is.null(tr)) {
    return(round(Q_loggumbel, 2))
  } else {
    match.tr <- match(tr, TR)
    if (is.na(match.tr)) {
      stop("TR not identified")
    } else {
      return(round(Q_loggumbel[match.tr], 2))
    }
  }
}
