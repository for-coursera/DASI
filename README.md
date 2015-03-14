Data Analysis and Statistical Inference by Dr. Mine Çetinkaya-Rundel
====================================================================


## Week 1

### Variables

  * Variables, numerical or categorical, quantitative or qualitative:
    * If the variable is numerical, further classify it as **continuous** or **discrete** based on whether or not the variable can take on an infinite number of values or only non-negative whole numbers, respectively.
    * If the variable is categorical, determine if it is **ordinal** based on whether or not the levels have a natural ordering.
  * Variables, associated or independent:
    *  Define **associated** variables as variables that show some relationship with one another. Further categorize this relationship as positive or negative association, when possible.
    * Define variables that are not associated as independent.
  * Variables, explanatory or response:
    * Identify the **explanatory** variable in a pair of variables as the variable suspected of affecting the other. The affected variable is called a **response** variable. However, note that labeling variables as explanatory and response does not guarantee that the relationship between the two is actually causal, even if there is an association identified between the two variables.

### Studies

  * Studies, observational or experimental
    * **Observational**: collect data in a way that doesn't directly interfere with how the data arise; such studies only _observe_ the data, only establish an association. Can be:
        * _retrospective_: uses past data
        * _prospective_: collects the data throughout the study
    * **Experimental**: randomly assign subject to treatments and establish casual connections; these studies lead _experiments_. These are the four principles of experimental design: 
        * control any possible confounders
        * randomize into treatment and control groups
        * replicate by using a sufficiently large sample or repeating the experiment
        * block any variables that might influence the response:
            * E.g., design an investigation whether energy gels help you run faster (treatment: energy gel, control: no energy gel)
                * however, energy gels might affect pro and amateur athletes differently. and thus, we need to block for pro status:
                    * divide the sample to pro and amateur
                    * randomly assign pro and amateur athletes to treatment and control groups
                    * pro and amateur athletes are equally represented in both groups

  * Experimental studies make use of:
    * **Random sampling**
        * Which can be of these properties:
            * **Simple random sampling**: Each subject in the population is equally likely to be selected. 
            * **Stratified sampling**: First divide the population into _homogenous_ strata (subjects within each stratum are similar, across strata are different), then randomly sample from within each strata (E.g., divide men and women into two separate stratas and then randomly sampling from these stratas). 
            * **Cluster sampling**: First divide the population into clusters (subjects within each cluster are _non-homogenous_, but clusters are similar to each other), then randomly sample a few clusters, and then randomly sample from within each cluster.
        * Or, can face the following bias:
            * **Convenience sample**: If individuals, who are easily accessible are more likely to be included in the sample
            * **Non-response**: If only a [non-random] fraction of the randomly sampled people respond to a survey such that sample is no longer representative of the popuplation
            * **Volunatry response**: If the sample consists of people, who volunteer to respond because they have strong oppinions on the issue.
    * **Random assignment** (or random placement) 
        * Which is an experimental technique for assigning subjects to different groups in an experiment (e.g., a treatment group versus a control group) using a chance procedure.
    * **Blinding**
        * Which can be
            * **Single**: If a subject is blinded
            * **Double**: If both tester and subject are blinded

  * If there sources of bias in a given study, these variables are called **confounding**.

### Results
  * Results, can be generalized to the population or not, and whether the results suggest correlation or causation between the quantities studied:
    * If random sampling has been employed in data collection, the results should be generalizable to the target population.
    * If random assignment has been employed in study design, the results suggest causality. 

### Analysis
  * Numerical variables
    * Pay attention to its shape, center, and spread, as well as any unusual observations
    * Note that there are three commonly used measures of center and spread: 
        * Center: **mean** (the arithmetic average), **median** (the midpoint), **mode** (the most frequent observation)
            * If median is larger than mean, then left skewness is possible
            * If mean is larger than median, then right skewness is possible
        * Spread: **standard deviation** (variability around the mean), **range** (max-min), **interquartile range** (middle 50% of the distribution)
    *  Identify the shape of a distribution as _symmetric_, _right skewed_, or _left skewed_, and _unimodal_, _bimodoal_, _multimodal_, or _uniform_.
    * Use _histograms_ and _box plots_ to visualize the shape, center, and spread of numerical distributions, and intensity maps for visualizing the spatial distribution of the data
    * Define a **robust statistic** (e.g. median, IQR) as a statistic that is not heavily affected by skewness and extreme outliers, and determine when such statistics are more appropriate measures of center and spread compared to other similar statistics
    * Recognize when transformations (e.g. log) can make the distribution of data more symmetric, and hence easier to model
    * Use _scatterplots_ for describing the relationship between two numerical variables, making sure to note the direction (positive or negative), form (linear or non-linear), and the strength of the relationship as well as any unusual observations that stand out (generally, we place the explanatory variable on the x-axis, and the response variable on the y-axis).
    * CORRELATION DOESN'T IMPLY CAUSATION!

  * Categorical variables
    * Use _frequency tables_ and bar plots to describe the distribution of one categorical variable
    * Use _contingency tables_ and segmented bar plots or mosaic plots to assess the relationship between two categorical variables
    * Use _side-by-side box plots_ for assessing the relationship between a numerical and a categorical variable

  * Note that an observed difference in sample statistics suggesting dependence between variables may be due to random chance, and that we need to use hypothesis testing to determine if this observed difference is too large to be attributed to random chance. Therefore, it's essential to set up _null_ and _alternative_ hypotheses for testing for independence between variables, and evaluate the data’s support for these hypotheses using a simulation technique:
    * Set a null and an alternative hypothesis
    * Simulate the experiment assuming that the null hypothesis is true
    * Thus, build a distribution of outcomes
    * Evaluated the probability (**p-value**) of observing an outcome at least as extreme as the one observed in the original data
    * And if this probability is sufficienty low, reject the null hypothesis in favor of the alternative


## Week 2

### Outcomes and events

  * Define the **experiment** as the situation involving probability that leads to results called **outcomes**.
  * Define the **outcome** as the result of a single trial of an experiment
  * Define the **event** as one or more outcomes of an experiment
  * Define the **probability** of an outcome as the proportion of times the outcome would occur if we observed the random process that gives rise to it an infinite number of times

  * Define **disjoint** (or **mutually exclusive**) **events** as events that cannot both happen at the same time: if A and B are disjoint, then `P(A and B) = 0`
  * Distinguish between **disjoint** and **independent** **events**:
    * If A and B are _independent_, then having information on A does not tell us anything about B (and vice versa)
    * If A and B are _disjoint_, then knowing that A occurs tells us that B cannot occur (and vice versa)
    * _Disjoint_ (mutually exclusive) events are always _dependent_ since if one event occurs we know the other one cannot

  * Define **complementary** outcomes as mutually exclusive outcomes of the same random process whose probabilities add up to 1: if A and B are complementary, then `P(A) + P(B) = 1`

  * Distinguish between **union of events (A or B)** and **intersection of events (A and B)**
  * Calculate the probability of _union of events_ using the (general) addition rule: 
    * If A and B _are not mutually exclusive_, then `P(A or B) = P(A) + P(B) − P(A and B)`
    * If A and B _are mutually exclusive_, `P (A or B) = P (A) + P (B)` (since for mutually exclusive events `P(A and B) = 0`)
  * Calculate the probability of _intersection of events_ using the multiplication rule: 
    * If A and B are _independent_, then `P(A and B) = P(A) × P(B)`
    * If A and B are _dependent_, then `P(A and B) = P(A|B) × P(B)`

### Probabilities of outcomes and events

  * Distinguish between **marginal**, **joint** and **conditional** probabilities
    * _Marginal_ probability: the probability of an event occurring (p(A)), it may be thought of as an unconditional probability.  It is not conditioned on another event. Example: the probability that a card drawn is red (p(red) = 0.5). Another example: the probability that a card drawn is a 4 (p(four)=1/13)
    * _Joint_ probability:  p(A and B).  The probability of event A and event B occurring.  It is the probability of the intersection of two or more events.  The probability of the intersection of A and B may be written p(A ∩ B). Example:  the probability that a card is a four and red =p(four and red) = 2/52=1/26.  (There are two red fours in a deck of 52, the 4 of hearts and the 4 of diamonds)
    * _Conditional_ probability:  p(A|B) is the probability of event A occurring, given that event B occurs. Example:  given that you drew a red card, what’s the probability that it’s a four (p(four|red))=2/26=1/13.  So out of the 26 red cards (given a red card), there are two fours so 2/26=1/13
  * Understand Bayes' Theorem: `P(A|B) = P(A and B) / P(B)` (thus, `P(A and B) = P(A|B) * P(B)` (see above))
    * [Breast cancer](http://sites.nicholas.duke.edu/statsreview/probability/jmc/) [test](http://www.bbc.com/news/magazine-28166019) example

### Normal distribution

  * Define the **standardized (Z) score** of a data point as the number of SDs it is away from the mean: `Z = (X − μ) / σ`
    * Obviously, Z-score of the _mean_ is always 0: `Z(μ) = 0`
  * Use the **Z score**
    * If the distribution is normal, to determine the _percentile_ score of a data point: 
        * Use `pnorm(X, μ, σ)` to detect the percentile (i.e., the probability that any event from this outcome lies below event X)
        * Use `qnorm(P, μ, σ)` to detect the value of event X, which is the upper border for all the events from this outcome, which occur with the probability P
        * Note the `lower.tail` option
    * Regardless of the shape of the distribution, to assess whether or not the particular observation is considered to be unusual (more than 2 standard deviations away from the mean)
    * Depending on the shape of the distribution determine whether the _median_ would have a negative (right skewed), positive (left skewed), or 0 (symmetrical) Z score keeping in mind that the mean always has a Z score of 0

### Binomial distribution

  * Determine if a random variable is **binomial** using the four conditions.
    * The trials are independent. 
    * The _number of trials_, `n`, is fixed. 
    * Each trial outcome can be classified as a success or failure. 
    * The _probability of a success_, `p`, is the same for each trial. 
  * Build a binomial distribution: `rbinom(n, SIZE, p)` (compare to `rnorm(SIZE, μ, σ)`)
  * Calculate the number of possible scenarios for obtaining `k` successes in `n`: `choose(n, k)`
  * Calculate the probability of a given number of successes in a given number of trials: `sum(dbinom(k:n, n, p))`
  * Calculate the expected number of successes in a given number of binomial trials `μ = n * p` and its standard deviation `σ = sqrt(n * p * (1 − p))` (see `bitono` function)
  * With a sufficiently large number of trials (n * p ≥ 10 and n * (1 − p) ≥ 10), use the normal approximation to calculate binomial probabilities (converting the binomial distribution in question to the normal distribution via `bitono` function)
