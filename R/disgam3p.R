#' Three-parameter Gamma distribution function for 24-hour maximum precipitation
#'
#' Three-parameter Gamma distribution function for 24-hour maximum precipitation for the study
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
#' @name disgam3p

disgam3p <- function(x, tr = NULL){
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

  x <- na.omit(x)
  x <- x * 1.13
  xm <- mean(x)
  dst <- sd(x)
  # cs <- e1071::skewness(x, type = 1) ### obs

  cs <- sum((x - xm)^3) / length(x) / (dst^3)

  coefs <- list(c = c(2.515517, 0.802853, 0.010328),
                d = c(1.432788, 0.189269, 0.001308))

  get_k_values <- function(n) {
    S <- seq(1, n)
    Px <- S / (n + 1)
    TM <- 1 / (1 - Px)
    Px_adj <- ifelse(Px >= 0.5, 1 - Px, Px)
    W <- sqrt(log(1 / Px_adj^2))
    k <- W - (coefs$c[1] + coefs$c[2] * W + coefs$c[3] * W^2) /
      (1 + coefs$d[1] * W + coefs$d[2] * W^2 + coefs$d[3] * W^3)
    list(k = k, Px = Px, TM = TM)
  }

  k_vals <- list(
    `2`   = get_k_values(99),    # TR = 2–100
    `3`   = get_k_values(199),   # TR = 200
    `4`   = get_k_values(499),   # TR = 500
    `5`   = get_k_values(999)    # TR = 1000
  )

  TR <- c(2, 5, 10, 20, 50, 100, 200, 500, 1000)
  vae <- c(k_vals$`2`$k[50],  # TR = 2
           k_vals$`2`$k[80],  # TR = 5
           k_vals$`2`$k[90],  # TR = 10
           k_vals$`2`$k[95],  # TR = 20
           k_vals$`2`$k[98],  # TR = 50
           k_vals$`2`$k[99],  # TR = 100
           k_vals$`3`$k[199], # TR = 200
           k_vals$`4`$k[499], # TR = 500
           k_vals$`5`$k[999]) # TR = 1000

  N <- length(x)
  gc <- cs / sqrt(N * (N - 1)) / (N - 2) * (1 + 8.5 / N)
  be <- (2 / gc)^2
  al <- dst / sqrt(be)
  y  <- xm - dst * sqrt(be)

  Q_pt <- al * be * (1 - 1 / (9 * be) + vae * sqrt(1 / (9 * be)))^3 + y

  if (is.null(tr)) {
    return(Q_pt)
  } else {
    match.tr <- match(tr, TR)
    if (is.na(match.tr)){
      stop("TR no identificado")
    } else {
      return(Q_pt[match.tr])
    }
  }
}
