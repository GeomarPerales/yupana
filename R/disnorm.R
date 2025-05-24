#' Normal distribution function for 24-hour maximum precipitation
#'
#' Normal distribution for 24-hour maximum precipitation for the study
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
#' @name disnorm

disnorm <- function(x, tr = NULL){
  if (is.data.frame(x)) {
    if (ncol(x) == 1) {
      x <- as.numeric(x)
      stop("Si 'x' es un data.frame, debe tener al menos dos columnas.")
    }
    x <- x[,2]

  } else if (is.numeric(x)) {
    x <- x
  } else {
    stop("El argumento 'x' debe ser un data.frame o un vector numerico.")
  }

  x <- as.numeric(x)
  x <- na.omit(x)

  if (length(x) == 0) {
    stop("No hay datos validos para procesar despues de eliminar NA.")
  }

  x <- x * 1.13
  xm <- mean(x)
  dst <- sd(x)
  cv <- dst / xm

  # Coeficientes de Wilson-Hilferty
  co <- 2.515517
  c1 <- 0.802853
  c2 <- 0.010328
  d1 <- 1.432788
  d2 <- 0.189269
  d3 <- 0.001308


  TR <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)
  Px <- 1 - 1/TR
  Px_adj <- ifelse(Px >= 0.5, 1 - Px, Px)
  W <- sqrt(log(1 / Px_adj^2))
  k <- W - (co + c1*W + c2*W^2) / (1 + d1*W + d2*W^2 + d3*W^3)

  Q_n <- xm + k * dst

  if (is.null(tr)) {
    return(Q_n)
  } else {
    match.tr <- match(tr, TR)
    if (is.na(match.tr)) {
      stop("TR no identificado")
    } else {
      return(Q_n[match.tr])
    }
  }
}
