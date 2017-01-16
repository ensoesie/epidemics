#' Fit Serfing's regression model with robust regression.
#'
#' Takes a time series and fits a robust regression model using the method described by Serfling (1963) in 'Methods for Current Statistical Analysis of Excess Pneumonia-Influenza Deaths'.
#' Helper function 'serf.design' constructs the design matrix according to parameters
#' 'trend.poly', 'harmonic', and 'period'.
#' See serf.fit function for example on how to fit standard model.
#'
#'
#' @param t1 is the integer index of the parameter 'ts' where model fitting begins; defaults to 1
#' @param t2 is the integer index of the parameter 'ts' where model fitting ends; defaults to length(ts)
#' @param t3 is the integer index of the time series where fitted values begin; defaults to min(t2+1,length(ts))
#' @param t4 is the integer index of the time series where fitted values end; defaults to min(t2+1,length(ts))
#' @param trend.poly is specified as an integer vector representing the degrees of the polynomial terms to be included; defaults to 1 i.e. single linear trend term
#' @param harmonic is also specified as an integer vector representing fractional multiples of the period e.g. 'harmonic=c(1,2)' specifies two harmonic pairs (two sine terms and two cosine terms) of one period and one half-period respectively
#' @param covar adds additional covariance factors but this feature not currently implemented; defaults to NULL
#' @param thresh is a numeric value specifying the epidemic threshold as a multiple of the residual standard deviation; defaults to 2.5
#' @param max.iter numeric value indicating the maximum number of iterations; default set to 10
#' @param period numeric value; defaults to 52.2
#' @param ylim range of y-axis for plotting
#' @param plot logical, generates plot of times series with threshold and epidemic observations; defaults to TRUE
#' @param output logical, outputs fitted values and epidemic observations; default set to TRUE
#' @param ts vector; response vector of numeric values representing the time series to be modeled. By default this is set to 'g.ts' the global time series which MUST be specified by the user before calling this function
#'
#'
#'@export
#'
#'

serf.robust <- function(t1 = 1, t2 = length(ts), t3 = min(t2 + 1, length(ts)), t4 = min(t2 + 1, length(ts)), trend.poly = 1,
    harmonic = c(1, 2, 3), covar = NULL, method = "a", c = 2.1, thresh = 2.1, max.iter = 10, period = g.period, ylim = range(ts[t1:t4]),
    plot = TRUE, plotagg = 0, main = g.main, output = TRUE, ts = g.ts) {
    ### Variable set-up

    N <- t2 - t1 + 1
    N.future <- t4 - t3 + 1
    k1 <- length(trend.poly)
    k2 <- length(trend.poly) + length(harmonic)
    y <- ts[t1:t2]
    y2 <- ts[t3:t4]

    ### 'serf.dat' and 'future.dat' contain trend and harmonic terms as specified by 'trend.poly' and 'harmonic'
    ### parameters, respectively.  Predictors for the model-fitting portion are contained in serf.dat, and predictors for
    ### the forecasting are contained in future.dat.  Both dataframes are built using auxiliary function serf.design()

    serf.dat <- serf.design(t1 = t1, t2 = t2, trend.poly = trend.poly, harmonic = harmonic, covar = covar, ts = y)
    future.dat <- serf.design(t1 = t3, t2 = t4, trend.poly = trend.poly, harmonic = harmonic, covar = covar, ts = y)

    ### Robust model fit for baseline

    base.fit <- lm(formula(data.frame(y, serf.dat)), data = serf.dat)
    w <- 1/(abs(base.fit$residuals))

    for (i in 1:max.iter) {
        next.fit <- lm(formula(data.frame(y, serf.dat)), data = serf.dat, weights = w)
        b <- as.numeric(next.fit$coef)
        res <- next.fit$residuals
        s <- median(abs(res)/0.6745)
        z <- res/s
        if (method == "a")
            w <- sin(z/(c/pi))/(z/(c/pi)) * (abs(z) < (c)) else w <- as.numeric(abs(z) < c)
    }

    fit <- next.fit
    final.w <- w
    rm(base.fit, next.fit, w)

    fitvals <- as.numeric(predict(fit, newdata = data.frame(y, serf.dat)))
    predvals <- as.numeric(predict(fit, newdata = data.frame(y2, future.dat)))
    sd.resid <- sd(subset(y - fitvals, final.w != 0))
    epi.thresh <- as.numeric((y - fitvals) > thresh * sd.resid)

    ### Try to decompose into trend and seasonal components

    if (k1 != 1) {
        trend <- rep(NA, N)
        seas <- rep(NA, N)
    } else {
        trend <- as.vector(rep(coef(fit)[1], N) + coef(fit)[2] * serf.dat[, 1])
        seas <- as.vector(as.matrix(serf.dat[, -1]) %*% coef(fit)[-c(1, 2)])
    }

    ### Plot, output, clean-up

    if (plot) {
        if (plotagg == 0) {
            tsplot(t1 = t1, t2 = t4, main = main, ylim = ylim, ts = ts)
            lines(1:N, fitvals, col = "blue", lwd = 3)
            lines(fitvals + thresh * sd.resid, lwd = 3, lty = 2)
            points(which(epi.thresh == 1), y[epi.thresh == 1], pch = 19, cex = 1.2, col = "red")
            lines((N + 1):(N + N.future), predvals, col = "red", lwd = 3)
        } else {
            pts <- agg.ts(y[1:(7 * (length(y)%/%plotagg))], plotagg)
            pfv <- agg.ts(fitvals[1:(7 * (length(y)%/%plotagg))], plotagg)
            pth <- agg.ts((fitvals + thresh * sd.resid)[1:(7 * (length(y)%/%plotagg))], plotagg)
            tsplot(main = main, ylim = range(pts), ts = pts)
            lines(pfv, col = "blue", lwd = 3)
            lines(pth, lwd = 3, lty = 2)
            points(which(pts > pth), pts[which(pts > pth)], pch = 19, cex = 1.2, col = "red")
        }
    }

    res <- list(c(t1, t2), as.numeric(coef(fit)), as.matrix(vcov(fit)), fitvals, (y - fitvals), trend, seas, epi.thresh,
        predvals, final.w, fit)
    names(res) <- c("times", "coef", "vcov", "fit.vals", "resids", "trend", "seas", "epi.thresh", "pred.vals", "final.w",
        "model.fit")

    if (output)
        print(res)
    serf.fit <- res
}

##############################
