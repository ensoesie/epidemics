#' Change time scale
#'
#' Takes a time series and collapses to different temporal units. See serf.fit function for example on how to fit model.
#'
#' @param ts response vector of numeric values representing the time series to be modeled. By default this is set to 'g.ts' the global time series which MUST be specified by the user before calling this function
#' @param len numeric value indicating vector length
#'
#' @export
#'
#'

agg.ts <- function(ts, len) {
    if ((length(ts)%%len) != 0)
        print(paste("WARNING: length of time series not divisible by ", len, "; last value may be truncated", sep = ""))
    n <- ceiling(length(ts)/len)

    res <- rep(0, n)
    for (i in 1:n) res[i] <- sum(ts[(len * (i - 1) + 1):(len * i)])
    return(res)
}
