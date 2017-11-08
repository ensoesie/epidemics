#'
#' Fit Serfing's regression model
#'
#' Takes a time series and fits a regression model using the method described by Serfling (1963) in 'Methods for Current Statistical Analysis of Excess Pneumonia-Influenza Deaths'.
#' Helper function 'serf.design' constructs the design matrix according to parameters
#' 'trend.poly', 'harmonic', and 'period'.
#'
#' @param ts response vector of numeric values representing the time series to be modeled. By default this is set to 'g.ts' the global time series which MUST be specified by the user before calling this function
#' @param t1 is the integer index of the parameter 'ts' where model fitting begins; defaults to 1
#' @param t2 is the integer index of the parameter 'ts' where model fitting ends; defaults to length(ts)
#' @param t3 is the integer index of the time series where fitted values begin; defaults to min(t2+1,length(ts))
#' @param t4 is the integer index of the time series where fitted values end; defaults to min(t2+1,length(ts))
#' @param trend.poly is specified as an integer vector representing the degrees of the polynomial terms to be included; defaults to 1 i.e. single linear trend term
#' @param harmonic is also specified as an integer vector representing fractional multiples of the period e.g. 'harmonic=c(1,2)' specifies two harmonic pairs (two sine terms and two cosine terms) of one period and one half-period respectively
#' @param period numeric value; defaults to 52.2
#' @param ylim range of y-axis for plotting
#' @param plot logical, generates plot of times series with threshold and epidemic observations; defaults to TRUE
#' @param output logical, outputs fitted values and epidemic observations; default set to TRUE
#' @param thresh is a numeric value specifying the epidemic threshold as a multiple of the residual standard deviation; defaults to 2.5
#'
#' @return epidemic numeric 0/1, 1 indicates epidemic period and 0 otherwise
#'
#'
#' @export
#'
#' @examples
#' library(splines)
#' library(MASS)
#' library(date)
#'
#' options(warn = 2)
#'
#' data(wa_who) # load WHO ILI data for Washington state
#' names(wa_who) <- c('id', 'ili', 'total', 'iliperc') # rename columns
#'
#' head(wa_who)
#'
#' wa_who$year <- floor(wa_who / 100)[1] # extract year
#' wa_who$week <- (wa_who %% 100)[1]   # extract week
#'
#' names(wa_who) <- c('id', 'ili', 'total', 'iliperc', 'year', 'week') # add column names
#' datelist <- sort(unique(wa_who$id)) # save sorted dates
#' year <- floor(datelist / 100)
#' week <- datelist %% 100
#'
#' ili <- rep(NA, length(datelist))
#' tot <- rep(NA, length(datelist))
#'
#' for (i in 1:length(datelist)) {
#'   ili[i] <- sum(wa_who$ili[wa_who$id == datelist[i]], na.rm = TRUE)
#'   tot[i] <- sum(wa_who$total[wa_who$id == datelist[i]], na.rm = TRUE)
#' }
#'
#' ili2.dat <- data.frame(id = (1:length(datelist)), year, week, ili = ili, total = tot)
#' ili2.dat$perc <- 100 * round(ili2.dat$ili / ili2.dat$total, 3)
#'
#' # Global data and plotting parameters
#' # Change as needed
#'
#' g.ts <- ili2.dat$perc
#' g.ts[is.nan(g.ts)] <- 0
#' g.period <- 52.2
#' g.mon <- 'Oct'
#' g.yr <- 1997 # data start year
#' g.xlab <- 'Calendar year'
#' g.ylab <- 'Percentage of samples positive'
#' g.main <- 'WHO'
#' g.dates <- TRUE
#'
#' resvec <- rep(0, length(g.ts))
#' stoptime <- length(g.ts) - 259
#' resmat <- matrix(NA, stoptime, length(g.ts))
#'
#' for (i in 1:stoptime) resmat[i, i:(i + 259)] <- serf.fit(g.ts, t1 = i, t2 = i + 259,
#'                                                          harmonic = 1, thresh = 1.96,
#'                                                          output = T)$epi.thresh
#'
#' for (i in 1:length(g.ts)) resvec[i] <- max(resmat[, i], na.rm = TRUE)
#'
#' ilires.dat <- data.frame(year = ili2.dat$year, week = ili2.dat$week,
#'                          perc = ili2.dat$perc, epidemic = resvec)
#' # write.csv(ilires.dat, ' ', row.names = FALSE) # write output to file
#'
#' # pdf('who-timeseries.pdf') # save plot to desired location
#' # tsplot(g.ts)
#' # points(which(resvec == 1), g.ts[resvec == 1], pch = 19, col = 'red', cex = 0.8)
#' # dev.off()

serf.fit <- function(ts = g.ts, t1 = 1, t2 = length(ts), t3 = min(t2 + 1, length(ts)), t4 = min(t2 + 1, length(ts)),
    trend.poly = 1, harmonic = c(1, 2), period = 52.2, ylim = range(ts[t1:t4]), plot = TRUE, output = TRUE, thresh = 2.5) {
    ### Variable set-up

    N <- t2 - t1 + 1
    N.future <- t4 - t3 + 1
    k1 <- length(trend.poly)
    k2 <- length(trend.poly) + length(harmonic)
    y <- ts[t1:t2]

    ### 'serf.dat' and 'future.dat' contain trend and harmonic terms as specified by 'trend.poly' and 'harmonic'
    ### parameters, respectively.  Predictors for the model-fitting portion are contained in serf.dat, and predictors for
    ### the forecasting are contained in future.dat.  Both dataframes are built using auxiliary function serf.design()

    serf.dat <- serf.design(t1 = t1, t2 = t2, trend.poly = trend.poly, harmonic = harmonic, ts = y)
    future.dat <- serf.design(t1 = t3, t2 = t4, trend.poly = trend.poly, harmonic = harmonic, ts = y)

    ### Two model fits: the first ('fit0') identifies potential epidemic observations; the second ('fit1') model the data,
    ### excluding previously identified epidemic periods.

    fit0 <- lm(formula(data.frame(y, serf.dat)), data = serf.dat)
    which.weeks <- as.logical(y < (fit0$fitted.values + thresh * sd(fit0$residuals)))

    ### Second model fit, excluding potential epidemic weeks

    fit1 <- lm(formula(data.frame(y, serf.dat)), data = serf.dat, subset = which.weeks)
    fitvals <- as.numeric(predict(fit1, newdata = serf.dat))
    predvals <- as.numeric(predict(fit1, newdata = future.dat))

    sd.resid <- sd(subset(y - fitvals, which.weeks))
    epi.thresh <- as.numeric((y - fitvals) > thresh * sd.resid)

    # print the baseline
    baseline <- thresh * sd.resid
    print(paste("Threshold", baseline), sep = "=")
    flush.console()

    ### Try to decompose into trend and seasonal components

    if (k1 != 1) {
        trend <- rep(NA, N)
        seas <- rep(NA, N)
    } else {
        trend <- as.vector(rep(coef(fit1)[1], N) + coef(fit1)[2] * serf.dat[, 1])
        seas <- as.vector(as.matrix(serf.dat[, -1]) %*% coef(fit1)[-c(1, 2)])
    }

    ### Plot, output, clean-up

    if (plot) {
        tsplot(main = "Original data with model fit", ylim = ylim, ts = y)
        lines(fitvals, col = "blue", lwd = 3)
        lines(fitvals + thresh * sd.resid, lwd = 3, lty = 2)
        points(which(epi.thresh == 1), y[epi.thresh == 1], pch = 19, cex = 1.2, col = "red")
    }

    res <- list(c(t1, t2), as.numeric(coef(fit1)), as.matrix(vcov(fit1)), fitvals, (y - fitvals), trend, seas, epi.thresh,
        predvals, fit1)
    names(res) <- c("times", "coef", "vcov", "fit.vals", "resids", "trend", "seas", "epi.thresh", "pred.vals", "model.fit")

    if (output)
        print(res)
    serf.fit <- res
}
