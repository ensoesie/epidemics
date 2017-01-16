#' Plot time series
#'
#'
#' tsplot takes a time series and produces a plot with x-axis measured
#' off in years.
#' Parameter 'ts' contains the time series of values to be plotted.
#'
#' Graphing parameters can be used as with the plot function.
#' Argument 'mon' specifies which month is written alongside the year on the
#' x-axis (defaults to 'Oct').  Argument 'yr' is the starting year (defaults to 1962).
#' Suppress dates from axis with argument 'dates=F'.
#'
#' @export
#'

tsplot <- function(ts, t1 = 1, t2 = length(ts), dates = TRUE, per = 52.2, mon = "Oct", yr = 1997, xlab = "Calendar year",
    ylab = "Percentage of samples positive", main = "WHO", lwd = 1, col = "black", xlim = (t1:t2), ylim = range(ts[t1:t2],
        na.rm = T)) {
    par(xaxt = "n")
    plot(1:(t2 - t1 + 1), ts[t1:t2], type = "l", xlab = xlab, ylab = ylab, main = main, lwd = lwd, col = col, ylim = ylim)
    par(xaxt = "s")
    num.yrs <- ceiling(length(g.ts)/per)
    labels <- NULL
    if (dates)
        for (i in 1:num.yrs) labels <- c(labels, paste(g.mon, g.yr + i - 1)) else labels <- rep("", num.yrs)
    axis(1, at = (0:(num.yrs - 1)) * per, labels = labels)
}
