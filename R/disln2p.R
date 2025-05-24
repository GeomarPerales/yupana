#' Two-parameter Log-Normal distribution function for 24-hour maximum precipitation
#'
#' Two-parameter Log-Normal distribution function for 24-hour maximum precipitation for the study
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
#' @name disln2p

disln2p <- function(x, tr = NULL){
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
  dst <- sd(x)
  cv <- dst / xm
  uy <- 0.5 * log(xm^2 / (1 + cv^2))
  dy <- sqrt(log(1 + cv^2))

  # Periodos de retorno
  TR <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)
  Px <- 1 - 1 / TR
  Px_adj <- ifelse(Px >= 0.5, 1 - Px, Px)

  # Coeficientes Wilson-Hilferty
  co <- 2.515517; c1 <- 0.802853; c2 <- 0.010328
  d1 <- 1.432788; d2 <- 0.189269; d3 <- 0.001308

  W <- sqrt(log(1 / Px_adj^2))
  k <- W - (co + c1*W + c2*W^2) / (1 + d1*W + d2*W^2 + d3*W^3)


  Q_ln2t <- exp(uy + k * dy)

  if (is.null(tr)) {
    return(Q_ln2t)
  } else {
    match.tr <- match(tr, TR)
    if (is.na(match.tr)) stop("TR no identificado")
    return(Q_ln2t[match.tr])
  }
}
