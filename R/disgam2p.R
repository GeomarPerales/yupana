#' Two-parameter Gamma distribution function for 24-hour maximum precipitation
#'
#' Two-parameter Gamma distribution function for 24-hour maximum precipitation for the study
#' of extreme precipitation events.
#'
#' @param x Vector or dataframe of daily maximum precipitation values
#' @param tr OPTIONAL, Vector of return period values, default parameter is NULL.
#'  Typical values for TR are 2, 5, 10, 20, 50, 100, 200, 500, 1000 years.
#'
#' @export
#'
#' @import stats
#' @importFrom stats na.omit sd qgamma
#'
#' @author Geomar Perales Apaico
#'
#' @name disgam2p

disgam2p <- function(x, tr = NULL){
  if(is.data.frame(x)){
    if(ncol(x) == 2){
      x <- as.numeric(x[,2])
    }else if(ncol(x) == 1){
      x <- as.numeric(x)
    } else {
      stop("values not defined")
    }
  } else if(is.numeric(x)){
    x <- x
  } else {
    stop("values not defined")
  }

  x <- as.numeric(x)
  x <- na.omit(x)

  x <- x * 1.13
  xm <- mean(x)
  s <- sd(x)

  co <- 2.515517
  c1 <- 0.802853
  c2 <- 0.010328
  d1 <- 1.432788
  d2 <- 0.189269
  d3 <- 0.001308

  TR <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)
  Px <- 1 - 1 / TR

  alpha <- (xm^2) / (s^2)
  beta <- (s^2) / xm

  # Escalado de los cuantiles (manual)
  Q_gamma <- qgamma(Px, shape = alpha, scale = beta)

  if (is.null(tr)) {
    return(Q_gamma)
  } else {
    match.tr <- match(tr, TR)
    if (is.na(match.tr)) {
      stop("TR not identified")
    } else {
      return(Q_gamma[match.tr])
    }
  }
}
