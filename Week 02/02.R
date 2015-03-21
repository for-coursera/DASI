# WEEK 02

SI.bitono <- function(n, p) {
# Convert the Binomial Distribution to the Normal Distribution
    mu <- n * p
    si <- sqrt(mu * (1 - p))
    vec <- c(mu, si)
    return(vec)
}
