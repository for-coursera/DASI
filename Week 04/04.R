# ANOVA

# Initial data ####
set.seed(23)
t1 <- rnorm(100)
set.seed(23)
t2 <- rnorm(76, 259, 8)
set.seed(23)
t3 <- rnorm(150, 800, 54)
set.seed(23)
t4 <- rnorm(123)
boxplot(t1, t2, t3, t4)
lst <- list(t1, t2, t3, t4)

# Step-by-step approach ####
# combine everything into a list
data <- list(t1, t2, t3)
# calculate the means for these vectors
mean_list <- lapply(data, mean); mean_list
# calculate the mean for all elements
mean_all <- mean(unlist(data)) ; mean_all
# calculate SST
st <- lapply(data, function(x) x <- (x - mean_all)^2) ; st
# sapply applies a function (sum) to a list (st) but returns its results as a vector, not as list (lapply)
sst <- sum(sapply(st, sum)) ; sst
# calculate SSG
sg <- lapply(mean_list, function(x) x <- (x - mean_all)^2) ; sg
counts <- lapply(data, length)
ssg <- sum(unlist(sg) * unlist(counts))
ssg
# calculate SSE
sse <- sst - ssg
# calculate degrees of freedom
dft <- length(unlist(data)) - 1
dfg <- length(data) - 1
dfe <- dft - dfg
# calculate means squares
msg <- ssg / dfg
mse <- sse / dfe
# calculate F statistic
f <- msg / mse
# compute p-value
pf(f, dfg, dfe, lower.tail = F)

# Factoring approach ####

# create a factor (categorical) variable
# more on those: http://www.ats.ucla.edu/stat/r/modules/factor_variables.htm
labels <- rep( c("t1","t2","t3"), c(length(t1), length(t2), length(t3)) ) ; labels
# create a unique variable with all the values
core <- c(t1, t2, t3) ; core
# combine these two
df <- data.frame(labels, core) ; df

# calculate means
# `aggregate` takes `core` data and separates it by levels 
# from the newly created list with `labels` factor variable
# `by` must be a list
mean_groups <- aggregate(df$core, by = list(df$labels), mean)$x ; mean_groups
mean_whole <- mean(core); mean_whole

# calculate SSG
length_groups <- aggregate( df$core, list(df$labels), length )$x
ssg <- sum( length_groups * (mean_groups - mean_whole)^2 ) ; ssg
# calculate SST
sst <- sum( (df$core - mean_whole)^2 ) ; sst
# calculate SSE
sse <- sst - ssg
# calculate degrees of freedom
dft <- length(df$labels) - 1
dfg <- length(levels(df$labels)) - 1
dfe <- dft - dfg
# calculate the required statistics
msg <- ssg / dfg
mse <- sse / dfe
f <- msg / mse
pf(f, dfg, dfe, lower.tail = F)

# Modeling approach ####
# create a factor (categorical) variable
labels <- rep( c("t1","t2","t3"), c(length(t1), length(t2), length(t3)) ) ; labels
# create a unique variable with all the values
core <- c(t1, t2, t3) ; core
# create a model
fit <- lm(core ~ labels)
# supply that model to `anova` function
anova(fit)

# Detecting differences ####

lst <- list(t1, t2, t3, t4)
significance_level <- 0.05
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
# sapply applies a function (sum) to a list (st) but returns its results as a vector, not as list (lapply)
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
