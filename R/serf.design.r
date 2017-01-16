#'
#' Helper function for Serfling models
#'
#' Given a time series and a specified number of harmonic and polynomial
#'  terms, this function builds the dataframe needed for the call to lm(). See serf.fit function for example on how to fit model.
#'
#' @param ts response vector of numeric values representing the time series to be modeled. By default this is set to 'g.ts' the global time series which MUST be specified by the user before calling this function
#' @param t1 is the integer index of the parameter 'ts' where model fitting begins; defaults to 1
#' @param t2 is the integer index of the parameter 'ts' where model fitting ends; defaults to length(ts)
#' @param trend.poly is specified as an integer vector representing the degrees of the polynomial terms to be included; defaults to 1 i.e. single linear trend term
#' @param harmonic is also specified as an integer vector representing fractional multiples of the period e.g. 'harmonic=c(1,2)' specifies two harmonic pairs (two sine terms and two cosine terms) of one period and one half-period respectively
#' @param period numeric value; defaults to 52.2
#' @param covar adds additional covariance factors but this feature not currently implemented; defaults to NULL
#'
#'
#' @export
#'
#' @return Data frame containing predictors for model fitting and forecasting

serf.design <- function(ts = g.ts, t1 = 1, t2 = length(ts), period = 52.2, trend.poly = 1, harmonic = c(1, 2), covar = NULL) {
    ### Variable set-up

    N <- t2 - t1 + 1
    k1 <- length(trend.poly)
    k2 <- k1 + 2 * length(harmonic)
    if (is.data.frame(covar))
        k3 <- k2 + dim(covar)[2] else k3 <- 0
    y <- as.numeric(ts)

    ### 'serf.dat' contains trend and harmonic terms as specified by 'trend.poly' and 'harmonic' parameters, respectively.

    serf.dat <- data.frame(matrix(0, nrow = N, ncol = k2))
    for (i in 1:k1) {
        serf.dat[, i] <- (t1:t2)^trend.poly[i]
        names(serf.dat)[i] <- paste("trend.poly.", i, sep = "")
    }
    for (i in 1:((k2 - k1)/2)) {
        temp <- k1 + (2 * i - 1)
        serf.dat[, temp] <- sin(2 * pi * (t1:t2)/(harmonic[i] * period))
        serf.dat[, temp + 1] <- cos(2 * pi * (t1:t2)/(harmonic[i] * period))
        names(serf.dat)[temp] <- paste("harmonic.sin.", i, sep = "")
        names(serf.dat)[temp + 1] <- paste("harmonic.cos.", i, sep = "")
        rm(temp)
    }
    if (k3 != 0) {
        serf.dat[, (k2 + 1):k3] <- covar[t1:t2, ]
        names(serf.dat[, (k2 + 1):k3]) <- names(covar)
    }

    return(serf.dat)
}
