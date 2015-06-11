# @knitr MAIN_MAIN

# SKEWED NORMAL DISTRIBUTION ####

# Generates a skewed normal distribution
SI.rnorm_skewed <- function(coef = 1, direction = "right", mu = 0, sigma = 1) {
    # generates a skewed normal distribution given the coef of skewness and its direction
    direction <- switch (direction,
                         right =  1,
                         left = -1)
    RNORM_SKEWED <- rnorm(5000, mu, sigma) + 3 * sigma
    RNORM_SKEWED <- RNORM_SKEWED[RNORM_SKEWED > 0 & RNORM_SKEWED < 6 * sigma]
    RNORM_SKEWED <- direction * (RNORM_SKEWED ^ (coef + 1))
    hist(RNORM_SKEWED)
    abline(v=mean(RNORM_SKEWED), col="red")
    abline(v=median(RNORM_SKEWED), col="blue")
    return(RNORM_SKEWED)
}

# BINOMIAL DISTRIBUTIONS ####

# Converts a Binomial Distribution to a Normal Distribution
SI.bitono <- function(n, p) {
    # converts a binomial distribution to a normal distribution 
    # given its success probability and a number of trials
    mu <- n * p
    si <- sqrt(mu * (1 - p))
    vec <- c(mu, si)
    return(vec)
}

# BAYESIAN TREE ####

# Calculates the Bayesian inference (1st level)
SI.bayes.tree <- function(outcomes) {
    compliments <- 1 - outcomes
    probs <- as.vector(rbind(outcomes, compliments))
    empty <- c("")
    probs <- c(empty, probs)
    rez <- c()
    rez[1] <- as.numeric(probs[2]) * as.numeric(probs[4])
    rez[2] <- as.numeric(probs[2]) * as.numeric(probs[5])
    rez[3] <- as.numeric(probs[3]) * as.numeric(probs[6])
    rez[4] <- as.numeric(probs[3]) * as.numeric(probs[7])
    return(rez)    
}

# STANDARD ERROR ####

# Caclulates the standard error for a normal distribution
SI.standart_error <- function(sigma, sample_size) {
    # caclulates the standard error given the variance and the size of a sample 
    tmp <- sigma^2/sample_size
    if (length(sigma) > 1) {
        se <- sqrt(sum(tmp))
    } else {
        se <- sqrt(tmp)
    }
    return(se)
}

# Caclulates the standard error for the distribution of proportions
SI.prop.standart_error <- function(population_proportion, sample_size, pool = FALSE) {
    # caclulates the standard error given the proportion of a population and the size of a sample 
    if (pool) {
        pooled_proportion <- sum(population_proportion) / sum(sample_size)
        return( sqrt( 
                (pooled_proportion * (1 - pooled_proportion) / sample_size[1]) + 
                (pooled_proportion * (1 - pooled_proportion) / sample_size[2])
                ) )
    }
    tmp <- (population_proportion * (1 - population_proportion)) / sample_size
    if (length(population_proportion) > 1) {
        se <- sqrt(sum(tmp))
    } else {
        se <- sqrt(tmp)
    }
    return(se)
}

# CONFIDENCE INTERVALS ####

# Important Reminder:
#
#   Commonly used confidence levels (CLs) are 90%, 95%, 98%, and 99%
#   "...With the confidence level A we state that the answer lies in B..."
#

# Calculates the CI for a population
SI.confidence_intreval <- function(confidence_level, mu, standard_error) {
    # calculates the CI of the population given the CL and the mean and the standard error of a sample
    value <- (1 - confidence_level) / 2
    return(c(qnorm(value, mu, standard_error), qnorm(confidence_level + value, mu, standard_error)))
}

# Calculates the needed sample size
SI.sample_size <- function(margin_of_error, sigma, confidence_level) {
    # calculates the needed sample size, given the standard deviation of a population, 
    # the desired margin of error and the confidence level of a projected sample
    z <- qnorm((1 - confidence_level) / 2)
    size <- ((z * sigma) / margin_of_error) ^ 2
    return(ceiling(size))
}

SI.prop.sample_size <- function(margin_of_error, population_proportion, confidence_level) {
    # calculates the needed sample size, given the standard deviation of a population, 
    # the desired margin of error and the confidence level of a projected sample    
    z <- qnorm((1 - confidence_level) / 2)
    return(ceiling(z^2 * population_proportion * (1 - population_proportion) / margin_of_error ^ 2))
}

# Important Reminder:
#
#   pnorm() takes values from a distribution as an input, and returns the probability with which that value
#           occurs in a given distribution
#   qnorm() - on the contrary - takes the probability as an input, and returns the value occurring with  
#           that probability in a given distribution
#

# HYPOTHESIS TESTING ####

# Calculates the z-score
SI.z.score <- function(test_value, mu, standard_error) {
    # calculates the z-score, given the proposed test value abd the mean and the standard error of a sample
    return((test_value - mu)/standard_error)
}

# Calculates the p-value
SI.pvalue <- function (null_hypothesis, standard_error, test_value_left = NULL, test_value_right = NULL) {
    # calculates the p-value, given the null hypothesis, the standard error and the test values
    if (is.null(test_value_left) && is.null(test_value_right)) {
        stop("No values to test provided")
    }
    if (!is.null(test_value_left) && is.null(test_value_right)) {
        return(pnorm((test_value_left - null_hypothesis)/standard_error))
    } 
    if (is.null(test_value_left) && !is.null(test_value_right)) {
        return(pnorm((test_value_right - null_hypothesis)/standard_error, lower.tail = F))
    } 
    if (!is.null(test_value_left) && !is.null(test_value_right)) {
        return( pnorm((test_value_left - null_hypothesis)/standard_error, lower.tail = T) 
                + pnorm((test_value_right - null_hypothesis)/standard_error, lower.tail = F) )
    } 
}

# Calculates p-value using standard deviation
SI.pvalue.sd <- function (null_hypothesis, sigma, sample_size, 
                          test_value_left = NULL, test_value_right = NULL) {
    # calculates p-value, given the mean and the variance of a population and the size of a sample
    if (is.null(test_value_left) && is.null(test_value_right)) {
        stop("No values to test provided")
    }
    if (!is.null(test_value_left) && is.null(test_value_right)) {
        return(pnorm((test_value_left - null_hypothesis)/(sigma/sqrt(sample_size))))
    }
    if (is.null(test_value_left) && !is.null(test_value_right)) {
        return(pnorm((test_value_right - null_hypothesis)/(sigma/sqrt(sample_size)), lower.tail = F))
    }
    if (!is.null(test_value_left) && !is.null(test_value_right)) {
        return( pnorm((test_value_left - null_hypothesis)/(sigma/sqrt(sample_size)), lower.tail = T)
                + pnorm((test_value_right - null_hypothesis)/(sigma/sqrt(sample_size)), lower.tail = F) )
    }
}

# BOOTSTRAPPING ####

# Generates a bootstrapping distribution
SI.boot <- function(data, boot_size, statistic) {
    # generates a bootstrapping distribution of a statistic given the original data and the required size
    boot <- rep(NA, boot_size)
    for (i in 1:boot_size) {
        boot[i] <- match.fun(statistic)(sample(data, length(data), replace = T))
    }
    return(boot)
}

# Calculates the CI for a bootstrapping distribution
SI.boot.confidence_interval <- function(data, confidence_interval = 0.9, method = "percentile") {
    # calculates the CI given the bootstrapping distribution and the required interval
    value <- (1 - confidence_interval) / 2
    if (method == "percentile") {
        quantile(data, c(value, 1-value))
    } else {
        if (method == "se") {
            mean <- mean(data)
            sd <- sd(data)
            v1 <- mean + qnorm(value) * sd
            v2 <- mean + qnorm(1-value) * sd
            rez <- c(v1, v2)
            names(rez) <- c(paste(value * 100, "%", sep=""), paste((1-value) * 100, "%", sep=""))
            return(rez)
        } else {
            stop("Invaid method")
        }
    }
}

# STUDENT'S T DISTRIBUTION ####

# Calculates the t-score
SI.t.score <- function(confidence_level = NULL, confidence_interval = NULL, sample_size) {
    # calculates the t-score, given the confidence level or the confidence interval & the size of a sample
    if (is.null(confidence_level) && is.null(confidence_interval)) {
        stop("Please, provide thr confidence level or the confidence interval")
    }
    if (is.null(confidence_interval)) {
        return( qt((1 - (1 - confidence_level) / 2), sample_size - 1) )
    }
    if (is.null(confidence_level)) {
        return( qt((1 - (1 - confidence_interval) / 2), sample_size - 1) )
    }
    
}

# Calculates the CI for a population
SI.t.confidence_intreval <- function(confidence_level, mean, standard_error, sample_size) {
    # calculates the CI of the population given the CL and the mean, the sd and the size of a sample
    tscore <- qt((1 - confidence_level) / 2, sample_size - 1)
    return( c(mean + tscore * standard_error, mean - tscore * standard_error ) )
}

# Calculates p-value
SI.t.pvalue <- function(null_hypothesis, standard_error, sample_size, 
                        test_value_left = NULL, test_value_right = NULL) {
    # calculates the p-value, given the null hypothesis, the standard error and the test values
    if (is.null(test_value_left) && is.null(test_value_right)) {
        stop("No values to test provided")
    }
    df <- min(sample_size - 1)
    if (!is.null(test_value_left) && is.null(test_value_right)) {
        return(pt((test_value_left - null_hypothesis)/standard_error, df))
    } 
    if (is.null(test_value_left) && !is.null(test_value_right)) {
        return(pt((test_value_right - null_hypothesis)/standard_error, df, lower.tail = F))
    } 
    if (!is.null(test_value_left) && !is.null(test_value_right)) {
        return( pt((test_value_left - null_hypothesis)/standard_error, df, lower.tail = T) 
                + pt((test_value_right - null_hypothesis)/standard_error, df, lower.tail = F) )
    } 
}

# ANOVA ####

# Compute analysis of variance (or deviance) tables for one or more fitted model objects
SI.anova <- function(list_of_distributions) {
    # compute analysis of variance tables given a list of distributions
    
    # analyze the data
    lst <- list_of_distributions
    len <- length(lst)
    dat <- unlist(lst)
    lab <- c()
    # transform the data according to the analysis
    for (i in 1:len) {
        vec <- rep(labels(lst)[i], length(lst[[i]]))
        lab <- append(lab, vec)
    }
    # create a model
    fit <- lm(dat ~ lab)
    # supply the model to `anova` function
    return(anova(fit))
}

SI.anova_pairwise <- function(list_of_distributions, significance_level = 0.05) {
    lst <- list_of_distributions
    len <- length(lst)
    # number of comparisons: (number of groups) * (number of groups - 1) / 2
    num <- len * (len - 1) / 2
    # correction Bonferroni: (significance level) / (number of comparisons)
    BC <- significance_level / num
    # to calculate SE for pairwise comparisons we need to calculate MSE:
    mean_list <- lapply(lst, mean); mean_list
    # calculate the mean for all elements
    mean_all <- mean(unlist(lst)) ; mean_all
    # calculate SST
    st <- lapply(lst, function(x) x <- (x - mean_all)^2)
    # sapply applies a function (sum) to a list (st), returns its results as a vector, not as list (lapply)
    sst <- sum(sapply(st, sum)) ; sst
    # calculate SSG
    sg <- lapply(mean_list, function(x) x <- (x - mean_all)^2) 
    counts <- lapply(lst, length)
    ssg <- sum(unlist(sg) * unlist(counts))
    ssg
    # calculate SSE
    sse <- sst - ssg
    # calculate degrees of freedom
    dft <- length(unlist(lst)) - 1
    dfg <- length(lst) - 1
    dfe <- dft - dfg
    # calculate means squares
    msg <- ssg / dfg
    mse <- sse / dfe ; mse
    # get names for every pair of distributions
    names <- combn(c(1:len), 2, simplify = F) ; names
    names <- as.character(names) ; names
    names <- gsub(":", ", ", gsub("c|\\(|\\)", "", names)) ; names
    # get list of lengths for every pair of distributions
    len_pairs <- combn(sapply(lst, length), 2, simplify = F)
    # separate it
    len_pairs1 <- lapply(len_pairs, "[[", 1)
    len_pairs2 <- lapply(len_pairs, "[[", 2)
    # produce a vector of SEs for every pair of distributions
    ses <- sqrt(sapply(len_pairs1, function(x) mse/x) + sapply(len_pairs2, function(x) mse/x))
    # produce a vector of the corresponding mean differences
    mean_pairs <- combn(sapply(lst, mean), 2, simplify = F)
    mean_diffs <- sapply(lapply(mean_pairs, rev), diff)
    # produce a vector of T statistics
    t <- mean_diffs / ses ; 
    # convert negative values to positive in place them to the right tail
    t <- abs(t)
    # calculate p-values for these pairs (for degrees of freedom we use dfe)
    p <- 2 * pt(t, dfe, lower.tail = F)
    # assign corresponding names
    names(p) <- names ; p
    # compare to Bonferroni correction
    fail <- p[p > BC]
    succ <- p[p < BC]
    print("The following list of distributions was provided")
    print(str(lst))
    if (!is.null(succ)) {
        print("Null hyphothesis successfully rejected for the pair(s): ")
        print(succ)    
    }
    if (!is.null(succ)) {
        print("Failed to reject null hyphothesis for the pair(s): ")
        print(fail)
    }
}

# CHI-SQUARE ####

# Calculates the chi-square
SI.chisq <- function(observed, expected) {
    # calculates the chi-square given the vectors of observed and expected outcomes
    if (length(observed) != length(expected)) {
        stop("Check data, please")
    }
    return(sum( (observed - expected)^2 / expected ))
}

# Calculates the p-value for a chi-square distribution
SI.chisq.pvalue <- function (chisq, number_of_levels) {
    # calculates the p-value given the chi-square statistic and the number of categorical levels
    # goodness of fit: one categorical variable, more than two levels
    return(pchisq(chisq, number_of_levels - 1, lower.tail = F))
}

SI.chisq.independence_test <- function (chisq, number_of_levels) {
    df <- (number_of_levels[1] - 1) * (number_of_levels[2] - 1)
    return(pchisq(chisq, df, lower.tail = F))
}

# MULTIPLE LINEAR REGRESSION MODELS ####

# Fits separate models for separate categories of a categorical variable on a scatterplot
SI.model.cat_separated <- function(model, ...){
    # fits separate models for separate categories of a categorical variable on a scatterplot
    # the code provided by DASI Team
    if(class(model)!="lm"){
        warning("Model must be the output of the function lm()")
    }
    
    if(length(model$xlevels)!=1){
        warning("Model must contain exactly one categorical predictor")
    }
    
    if(length(model$coef)-length(model$xlevels[[1]])!=1){
        warning("Model must contain exactly one non-categorical predictor")
    }
    
    palette <- c("#E69F00", "#56B4E9", "#D55E00", "#009E73", "#CC79A7", "#F0E442", "#0072B2")
    
    baseIntercept <- model$coef[1]
    nLines <- length(model$xlevels[[1]])
    intercepts <- c(baseIntercept, rep(0, nLines-1))
    indicatorInd <- c(1, rep(0, nLines)) # used to find slope parameter by process of elimination
    
    for(i in 1:(nLines-1)){
        indicatorName <- paste(names(model$contrasts),model$xlevels[[1]][1+i], sep = "")
        intercepts[i+1] <- baseIntercept + model$coef[names(model$coef)==indicatorName]
        indicatorInd <- indicatorInd + (names(model$coef)==indicatorName)
    }
    
    slope <- model$coef[!indicatorInd]
    
    num_pred = which(names(model$model[,-1]) != names(model$xlevels)) + 1
    cat_pred = which(names(model$model[,-1]) == names(model$xlevels)) + 1
    
    model$model$COL = NA
    model$model$PCH = NA
    for(i in 1:nLines){
        model$model$COL[model$model[,cat_pred] == levels(model$model[,cat_pred])[i]] = adjustcolor(palette[i],0.40)
        model$model$PCH[model$model[,cat_pred] == levels(model$model[,cat_pred])[i]] = i+14
    }
    
    plot(model$model[,1] ~ jitter(model$model[,num_pred]), col = model$model$COL, pch = model$model$PCH,
         ylab = names(model$model)[1],
         xlab = names(model$model)[num_pred])
    
    for(j in 1:nLines){
        abline(intercepts[j], slope, col = palette[j], lwd = 2, ...)
    }
    
    if(slope > 0){legend_pos = "bottomright"}
    if(slope < 0){legend_pos = "topleft"}  
    
    legend(legend_pos, col = palette[1:nLines], lty = 1, legend = levels(model$model[,cat_pred]), lwd = 2)
}
