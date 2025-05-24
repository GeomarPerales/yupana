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

dislgum <- function(x, tr = NULL){
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
  xlog <- log(x)
  xm <- mean(xlog)
  dst <- sd(xlog)

  n <- length(xlog)
  Px1 <- c(0.5, 0.33, 0.25, 0.2, 0.167, 0.143, 0.125, 0.111, 0.1, 0.091)
  Px2 <- c(0.99, 0.961, 0.938, 0.922, 0.909, 0.899, 0.89, 0.882, 0.875, 0.868)
  Px3 <- c(0.995, 0.9756, 0.9608, 0.95, 0.9411, 0.9336, 0.9272, 0.9216, 0.9167, 0.9123)
  Px4 <- c(0.998, 0.9872, 0.9786, 0.9718, 0.9662, 0.9615, 0.9575, 0.954, 0.9508, 0.9479)
  Px5 <- c(0.999, 0.9936, 0.987, 0.9822, 0.9783, 0.9749, 0.9719, 0.9693, 0.967, 0.9649)

  TR <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11)

  Px <- switch(
    as.character(n),
    "10" = Px1,
    "20" = Px2,
    "30" = Px3,
    "40" = Px4,
    "50" = Px5
  )

  Y_gt <- -log(-log(Px))

  a1 <- 1.2825 / dst
  u <- xm - 0.45 * dst

  Q_log <- (Y_gt / a1) + u

  Q_gt <- exp(Q_log)

  if (is.null(tr)) {
    return(round(Q_gt[match(tr, TR)], 2))
  } else {
    match.tr <- match(tr, TR)
    if (is.na(match.tr)) {
      stop("TR not identified")
    } else {
      return(Q_gt[match.tr])
    }
  }

}
