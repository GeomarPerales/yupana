#' Three-parameter Log-Normal distribution function for 24-hour maximum precipitation
#'
#' Three-parameter Log-Normal distribution function for 24-hour maximum precipitation for the study
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
#' @name disln3p

disln3p <- function(x, tr = NULL){
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

  calc_weibull <- function(S, co, c1, c2, d1, d2, d3) {
    Px <- seq(1, length(S)) / (length(S) + 1)
    Px[Px >= 0.5] <- 1 - Px[Px >= 0.5]
    W <- sqrt(log(1 / Px^2))
    k <- W - (co + c1 * W + c2 * W^2) / (1 + d1 * W + d2 * W^2 + d3 * W^3)
    return(k)
  }

  sg <- sum((x - xm)^3) / length(x)

  # Coeficientes Weibull
  co <- 2.515517
  c1 <- 0.802853
  c2 <- 0.010328
  d1 <- 1.432788
  d2 <- 0.189269
  d3 <- 0.001308

  TR <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)

  sample_sizes <- c(59, 99, 199, 499, 999)

  k_values <- sapply(sample_sizes, function(n) calc_weibull(seq(1, n), co, c1, c2, d1, d2, d3))

  # Seleccionamos los valores correspondientes de k para los periodos de retorno
  vae <- c(k_values[[2]][50], k_values[[2]][80], k_values[[2]][90], k_values[[2]][95],
           k_values[[2]][98], k_values[[2]][99], k_values[[3]][199], k_values[[4]][499], k_values[[5]][999])

  g <- length(x)^2 * sg / (length(x) - 1) / (length(x) - 2) / dst^3
  cs <- g
  W <- (-g + sqrt(g^2 + 4)) * 0.5
  Z2 <- (1 - W^(2/3)) / W^(1/3)
  dy2 <- sqrt(log(Z2^2 + 1))
  uy2 <- log(dst / Z2) - 0.5 * log(Z2^2 + 1)
  xo <- xm - dst / Z2

  Q_ln3t <- sapply(vae, function(v) xo + exp(uy2 + v * dy2))

  if (is.null(tr)) {
    return(Q_ln3t)
  } else {
    match.tr <- match(tr, TR)
    if (is.na(match.tr)) {
      stop("TR not identified")
    } else {
      return(Q_ln3t[match.tr])
    }
  }
}
