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
        * block any variables that might influence the response

  * Experimental studies make use of:
    * **Random sampling**
        * Which can be of these properties:
            * **Simple random sampling**: Each subject in the population is equally likely to be selected. 
            * **Stratified sampling**: First divide the population into homogenous strata (subjects within each stratum are similar, across strata are different), then randomly sample from within each strata. 
            * **Cluster sampling**: First divide the population into clusters (subjects within each cluster are non-homogenous, but clusters are similar to each other), then randomly sample a few clusters, and then randomly sample from within each cluster.
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
        * Spread: **standard deviation** (variability around the mean), **range** (max-min), **interquartile range** (middle 50% of the distribution)
    *  Identify the shape of a distribution as _symmetric_, _right skewed_, or _left skewed_, and _unimodal_, _bimodoal_, _multimodal_, or _uniform_.
    * Use _histograms_ and _box plots_ to visualize the shape, center, and spread of numerical distributions, and intensity maps for visualizing the spatial distribution of the data
    * Define a **robust statistic** (e.g. median, IQR) as a statistic that is not heavily affected by skewness and extreme outliers, and determine when such statistics are more appropriate measures of center and spread compared to other similar statistics
    * Recognize when transformations (e.g. log) can make the distribution of data more symmetric, and hence easier to model
    * Use _scatterplots_ for describing the relationship between two numerical variables, making sure to note the direction (positive or negative), form (linear or non-linear), and the strength of the relationship as well as any unusual observations that stand out

  * Categorical variables
    * Use _frequency tables_ and bar plots to describe the distribution of one categorical variable
    * Use _contingency tables_ and segmented bar plots or mosaic plots to assess the relationship between two categorical variables
    * Use _side-by-side box plots_ for assessing the relationship between a numerical and a categorical variable

  * Note that an observed difference in sample statistics suggesting dependence between variables may be due to random chance, and that we need to use hypothesis testing to determine if this observed difference is too large to be attributed to random chance. Therefore, it's essential to set up _null_ and _alternative_ hypotheses for testing for independence between variables, and evaluate the data’s support for these hypotheses using a simulation technique.


## Week 2
