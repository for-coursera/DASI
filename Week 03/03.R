# WEEK 03

# CENTRAL LIMIT THEOREM ####

# Applying CLT 
# generate a population
pop <- rnorm(10000, 0, 20)
# compute its mean and sd
mu <- mean(pop) ; mu ; sd <- sd(pop) ; sd
# generate a list of samples of that population (replacement = FALSE, size < 10%)
sample_dst <- list() ; sample_size <- 300 ; for (i in 1:50) { sample_dst[[i]] <- sample(pop, sample_size) }
# compute a sampling distribution of means of that population
sampling_dst <- as.numeric(lapply(sample_dst, mean))
# compute the mean and the sd of the sampling distribution
mu_dst <- mean(sampling_dst) ; mu_dst ; sd_dst <- sd(sampling_dst) ; sd_dst
# compare statistics of the sampling distribution with the statistics of the population
abs(mu - mu_dst) ; abs((sd / sqrt(sample_size)) - sd_dst)
# build histograms
hist(pop) ; hist(sampling_dst)
 

# CONFIDENCE INTERVALS ####

# Important Reminder:
#
#   Commonly used confidence levels (CLs) are 90%, 95%, 98%, and 99%
#   "...With the confidence level A we state that the answer lies in B..."
#

# Calculates the CI of the population
SI.confidence_intreval <- function(confidence_level, sample_size, mu = 0, sigma = 1) {
# calculates the CI of the population given the CL and the sample's size, mean and variance
    value <- (1 - confidence_level) / 2
    standard_error <- sigma / sqrt(sample_size)
    return(c(qnorm(value, mu, standard_error), qnorm(confidence_level + value, mu, standard_error)))
}

# Calculates the needed sample size
SI.sample_size <- function(sigma, margin_of_error, confidence_level) {
# calculates the needed sample size, given the standard deviation of the population 
# and the desired margin of error and the confidence level of the projected sample
    z <- qnorm((1 - confidence_level) / 2)
    size <- ((z * sigma) / margin_of_error) ^ 2
    return(ceiling(size))
}

# Important Reminder:
#
#   pnorm() takes values from a distribution as an input, and returns the probability with which that value
#           occurs in a given distribution
#   qnorm() - on the contrary - takes the probability as an input, and returns the value occurring with  
#           that probability in a given distribution
#

# HYPOTHESIS TESTING ####

# Calculates z-score
SI.zscore <- function(mu, sigma, sample_size, test_value) {
# calculates z-score, given mean and variance of the population and mean and size of the sample
    return((test_value - mu)/(sigma/sqrt(sample_size)))
}

# Calculates p-value
SI.pvalue <- function (mu, sigma, sample_size, test_value_left = NULL, test_value_right = NULL) {
# calculates p-value, given mean and variance of the population and mean and size of the sample
    if (is.null(test_value_left) && is.null(test_value_right)) {
        stop("No values to test provided")
    }
    if (!is.null(test_value_left) && is.null(test_value_right)) {
        return(pnorm((test_value_left - mu)/(sigma/sqrt(sample_size))))
    } 
    if (is.null(test_value_left) && !is.null(test_value_right)) {
        return(pnorm((test_value_right - mu)/(sigma/sqrt(sample_size)), lower.tail = F))
    } 
    if (!is.null(test_value_left) && !is.null(test_value_right)) {
        return( pnorm((test_value_left - mu)/(sigma/sqrt(sample_size)), lower.tail = T) 
                   + pnorm((test_value_right - mu)/(sigma/sqrt(sample_size)), lower.tail = F) )
    } 
}

