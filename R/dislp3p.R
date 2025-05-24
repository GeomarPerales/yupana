#' Three-parameter Log-Pearson distribution function for 24-hour maximum precipitation
#'
#' Three-parameter Log-Pearson distribution function for 24-hour maximum precipitation for the study
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
#' @name dislp3p

dislp3p <- function(x, tr = NULL){
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


  co <- 2.515517
  c1 <- 0.802853
  c2 <- 0.010328
  d1 <- 1.432788
  d2 <- 0.189269
  d3 <- 0.001308


  calc_k <- function(n) {
    Px <- seq(1, n) / (n + 1)
    TM <- 1 / (1 - Px)
    Px_adj <- ifelse(Px >= 0.5, 1 - Px, Px)
    W <- sqrt(log(1 / Px_adj^2))
    k <- W - (co + c1 * W + c2 * W^2) / (1 + d1 * W + d2 * W^2 + d3 * W^3)
    list(k = k, Px = Px, TM = TM)
  }


  out2 <- calc_k(99)
  out3 <- calc_k(199)
  out4 <- calc_k(499)
  out5 <- calc_k(999)


  TR <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)
  vae <- c(out2$k[50], out2$k[80], out2$k[90], out2$k[95], out2$k[98], out2$k[99],
           out3$k[199], out4$k[499], out5$k[999])


  x <- x[,2]
  x <- as.numeric(x)
  x <- na.omit(x)
  x <- x * 1.13

  lnQ <- log(x)
  N <- length(x)
  xm_lp <- mean(lnQ)
  ds_lp <- sd(lnQ)
  cv_lp <- xm_lp / ds_lp


  sg_lp <- sum((lnQ - xm_lp)^3) / N
  g_lp <- N^2 * sg_lp / ((N - 1) * (N - 2) * ds_lp^3)
  cs_lp <- g_lp


  gc_lp <- cs_lp / (sqrt(N * (N - 1)) / (N - 2) * (1 + 8.5 / N))
  be_lp <- (2 / gc_lp)^2
  sc <- ds_lp * sqrt(N / (N - 1))
  al_lp <- sc / sqrt(be_lp)
  y_lp <- xm_lp - al_lp * be_lp


  Q_lpt <- exp(al_lp * be_lp * (1 - 1 / (9 * be_lp) + vae * sqrt(1 / (9 * be_lp)))^3 + y_lp)

  if (is.null(tr)) {
    return(Q_lpt)
  } else {
    match.tr <- match(tr, TR)
    if (is.na(match.tr)) {
      stop("TR not identified")
    } else {
      return(Q_lpt[match.tr])
    }
  }
}
