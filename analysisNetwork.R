# Load libraries.
library(haven)
library(dplyr)
library(tidyverse)
library(sjPlot)
library(lavaan)
library(powerly)
library(psych)
library(corrplot)
library(GGally)
library(semPlot)
library(bootnet)
library(qgraph)
library(NetworkComparisonTest)
library(BGGM)
library(CliquePercolation)
library(EGAnet)

# Number of cores available on the machine.
cores <- parallel::detectCores() - 1


# #region Sample size analysis.

# Set a seed.
set.seed(2001)

# Hypothesized true model.
model <- generate_model(
  type = "ggm",
  nodes = 17,
  density = 0.4
)

# Plot the true model.
qgraph::qgraph(model, layout = "spring")

#Run the power analysis.
results <- powerly(
    range_lower = 1000,
    range_upper = 1300,
    samples = 50,
    replications = 50,
    model = "ggm",
    model_matrix = model,
    measure = "sen",
    statistic = "power",
    measure_value = 0.8,
    statistic_value = 0.8,
    boots = 10000,
    cores = cores
)

# Run validation.
validation <- validate(
    method = results,
    replications = 3000,
    cores = cores
)

# Plot the analysis results.
plot(results)

# Plot the validation results.
plot(validation)

# #endregion


# #region Data pre-processing.

# Set the working directory accordingly.
setwd("./data")

# Load the data.
data <- read_sav("THORESCI_Selection_Vida_anon.sav")

# Glimpse the data.
data %>%
    glimpse()

# Visualize the labels.
data %>%
    view_df()

# Process the data.
data <- data %>%
    # Lowercase the variable names.
    rename_all(tolower) %>%
    # Rename the unhealthy lifestyle variables.
    rename(
        smoking = smoking_status,
        alcohol = alcohol_consumption_bl,
        bmi = bmi_bl,
        diet = diet_restrictfatssalt,
        activity = physical_activity_t0
    ) %>%
    # Recode variables, s.t., higher values indicate an unhealthy lifestyle.
    mutate(
        # Whether the participant restricts fats and salt in their diet.
        indulgence = 9 - diet,

        # Whether the participant is physically active.
        sedentary = 5 - activity,

        # Whether the participant currently smokes.
        smoking = factor(
            ifelse(smoking == 1, 1, 0),
            levels = c(0, 1),
            labels = c("no", "yes")
        ),

        # Whether the participant currently consumes alcohol.
        alcohol = factor(
            ifelse(alcohol == 1, 1, 0),
            levels = c(0, 1),
            labels = c("no", "yes")
        )
    ) %>%

    # Select relevant variables.
    dplyr::select(
        # Select the psychological variables.
        dplyr::starts_with(match = c("phq", "gad")),

        # Select the unhealthy lifestyle variables.
        smoking, alcohol, bmi, sedentary, indulgence
    ) %>%

    # Filter for complete cases on the unhealthy lifestyle variables.
    filter(
        complete.cases(smoking, alcohol, bmi, sedentary, indulgence)
    )

# Glimpse the processed data.
data %>%
    glimpse()

# Subset the data for the first time point.
data_t1 <- data %>%
    # Select the variables for the first time point.
    dplyr::select(ends_with(match = "t0")) %>%

    # Remove the time suffix.
    rename_all(~ gsub(pattern = "_t0", replacement = "", x = .))

# Subset the data for the second time point.
data_t2 <- data %>%
    # Select the variables for the second time point.
    dplyr::select(ends_with(match = "t2")) %>%

    # Remove the time suffix.
    rename_all(~ gsub(pattern = "_t2", replacement = "", x = .))

# Subset the data for the third time point.
data_t3 <- data %>%
    # Select the variables for the third time point.
    dplyr::select(ends_with(match = "t3")) %>%

    # Remove the time suffix.
    rename_all(~ gsub(pattern = "_t3", replacement = "", x = .))


# Subset the unhealthy lifestyle variables.
data_unhealthy <- data %>%
    # Select the unhealthy lifestyle variables.
    dplyr::select(smoking, alcohol, bmi, sedentary, indulgence)

# #endregion


# #region Descriptive statistics.

# Describe the data.
data %>%
    describe()

# Correlation plotting function to avoid code repetition.
plot_correlations <- function(data, time_point) {
    # Compute the correlations.
    cor(x = data, use = "complete.obs") %>%
    # Plot the correlations.
    corrplot(
        corr = .,
        method = "color",
        type = "lower",
        tl.col = "black",
        mar = c(0, 0, 1.5, 0),
        title = ("Time Point ")
    )
}

# Split the plot area.
layout(matrix(1:3, nrow = 1))

# Visualize the data for the first time point.
data_t1 %>%
    plot_correlations(time_point = 1)

# Visualize the data for the second time point.
data_t2 %>%
    plot_correlations(time_point = 2)

# Visualize the data for the third time point.
data_t3 %>%
    plot_correlations(time_point = 3)

# Reset the plot area.
layout(1)

# Visualize the unhealthy lifestyle variables.
data_unhealthy %>%
    ggpairs()

# #endregion


# #region Confirmatory factor analysis.

# Order the categorical variables for `lavaan`.
data_cfa <- data_unhealthy %>%
    # Declare the categorical variables as ordered.
    mutate(
        smoking = ordered(smoking),
        alcohol = ordered(alcohol)
    ) %>%

    # Obtain a standard data frame.
    as.data.frame()

# Specify `lavaan` measurement model for unhealthy lifestyle.
model <- "
    # Unhealthy lifestyle factor.
    unhealthy =~ smoking + alcohol + bmi + sedentary + indulgence
"

# Fit the model.
fit <- cfa(
    model = model,
    data = data_cfa,
    std.lv = TRUE,
    estimator = "WLSMV",
)

# Summarize the model.
summary(fit, standardized = TRUE, fit.measures = TRUE)

# Plot the model.
semPaths(fit, what = "path", whatLabels = "est", intercepts = FALSE)

# Obtain the factor scores.
unhealthy <- lavPredict(fit, type = "lv", transform = TRUE)

# Describe the factor scores.
describe(unhealthy)

# Add the unhealthy lifestyle factor to the first time point data.
data_t1 <- data_t1 %>%
    mutate(unhealthy = unhealthy)

# Add the unhealthy lifestyle factor to the second time point data.
data_t2 <- data_t2 %>%
    mutate(unhealthy = unhealthy)

# Add the unhealthy lifestyle factor to the third time point data.
data_t3 <- data_t3 %>%
    mutate(unhealthy = unhealthy)

# #endregion

## Region correlations between the factor and the surveys total score

#Calculate total depression survey score- PHQ9, first time point
PHQ_Total1 <- rowSums(data_t1[ , c(1:9)])
#Calculate Correlation time 1
cor.test(PHQ_Total1, unhealthy, na.rm=TRUE) 

#Calculate total anxiety survey score- GAD7, first time point
GAD_Total1 <- rowSums(data_t1[ , c(10:16)])
#Calculate Correlation time 1
cor.test(GAD_Total1, unhealthy, na.rm=TRUE)

## endregion

# #region Network analysis via regularization approach.

# Regularized network for the first time point.
fit_network_regularization_t1 <- estimateNetwork(
    data = data_t1,
    default = "EBICglasso",
    corMethod = "cor_auto"
)

# Regularized network for the second time point.
fit_network_regularization_t2 <- estimateNetwork(
    data = data_t2,
    default = "EBICglasso",
    corMethod = "cor_auto"
)

# Regularized network for the third time point.
fit_network_regularization_t3 <- estimateNetwork(
    data = data_t3,
    default = "EBICglasso",
    corMethod = "cor_auto"
)

# Specify the nodes to be included in the network.
nodes <- colnames(data_t1)[-ncol(data_t1)]
#nodes <- colnames(data_t1)
# Function for sub-setting nodes.
subset_nodes <- function(graph) {
    colnames(graph) %in% nodes
}
#Legend nodes
#Node names long
long_names <- c("Low interest", # PHQ1
                "Depressed Mood", # PHQ2
                "Disturbed Sleep", # PHQ3
                "Fatigue", # PHQ4
                "Disturbed Appetite", # PHQ5
                "Guilt", # PHQ6  
                "Concentration Problems", # PHQ7
                "Psychomotor Disturbances", # PHQ8 !!
                "Self-harm Thoughts", # PHQ9
                "Anxiety", # GAD1
                "Excessive Worry", # GAD2
                "Generalized Worry", # GAD3
                "Trouble Relaxing", # GAD4
                "Restlessness", # GAD5
                "Irritable Mood", # GAD6
                "Dread" # GAD7
)
#Items matched to constructs
item_groups_names <- list("Depression"=1:9,"Anxiety"=10:16)

# Extract the graph for the first time point.
graph_regularization_t1 <- fit_network_regularization_t1$graph %>%
    subset(subset_nodes(.), subset_nodes(.))

# Extract the graph for the second time point.
graph_regularization_t2 <- fit_network_regularization_t2$graph %>%
    subset(subset_nodes(.), subset_nodes(.))

# Extract the graph for the third time point.
graph_regularization_t3 <- fit_network_regularization_t3$graph %>%
    subset(subset_nodes(.), subset_nodes(.))

# Split the plot area.
layout(matrix(1:3, nrow = 1))

# Plot and store the regularized network for the first time point.
plot_network_regularization_t1 <- qgraph(
    input = graph_regularization_t1,
    labels = nodes,
    layout = "circle",
    theme = "colorblind",
    title = "Time Point 1",
    legend = FALSE,
    groups = item_groups_names
)

# Plot and store the regularized network for the second time point.
plot_network_regularization_t2 <- qgraph(
    input = graph_regularization_t2,
    labels = nodes,
    layout = "circle",
    theme = "colorblind",
    title = "Time Point 2",
    legend = FALSE,
    groups = item_groups_names
)

# Plot the regularized network for the third time point.
plot_network_regularization_t3 <- qgraph(
    input = graph_regularization_t3,
    labels = nodes,
    nodeNames = long_names,
    layout = "circle",
    theme = "colorblind",
    title = "Time Point 3",
    groups = item_groups_names,
    legend = FALSE,
    legend.cex = .55
)

# Reset the plot area.
layout(1)

# #endregion


# #region Network analysis via Bayesian approach.

# Function to estimate networks via `BGGM` and avoid repetition.
estimate_network_bayesian <- function(data) {
    # Add one to the ordinal variables to shift the scale and avoid zeros.
    mutate(.data = data, across(!unhealthy, ~ .x + 1)) %>%

    # Convert to regular data frame.
    as.data.frame() %>%

    # Estimate the network.
    estimate(
        Y = .,
        formula = ~ unhealthy,
        type = "ordinal",
        impute = FALSE
    )
}

# Bayesian network for the first time point.
fit_network_bayesian_t1 <- data_t1 %>%
    estimate_network_bayesian()

# Bayesian network for the second time point.
fit_network_bayesian_t2 <- data_t2 %>%
    estimate_network_bayesian()

# Bayesian network for the third time point.
fit_network_bayesian_t3 <- data_t3 %>%
    estimate_network_bayesian()


# Select the graph for the first time point.
graph_bayesian_t1 <- BGGM::select(fit_network_bayesian_t1)

# Select the graph for the second time point.
graph_bayesian_t2 <- BGGM::select(fit_network_bayesian_t2)

# Select the graph for the third time point.
graph_bayesian_t3 <- BGGM::select(fit_network_bayesian_t3)


# Split the plot area.
layout(matrix(1:3, nrow = 1))

# Plot the Bayesian graph for the first time point.
plot_network_bayesian_t1 <- qgraph(
    input = graph_bayesian_t1$pcor_adj,
    labels = nodes,
    layout = "circle",
    theme = "colorblind",
    title = "Time Point 1",
    groups = item_groups_names
)

# Plot the Bayesian graph for the second time point.
plot_network_bayesian_t2 <- qgraph(
    input = graph_bayesian_t2$pcor_adj,
    labels = nodes,
    layout = "circle",
    theme = "colorblind",
    title = "Time Point 2",
    groups = item_groups_names
)

# Plot the Bayesian graph for the third time point.
plot_network_bayesian_t3 <- qgraph(
    input = graph_bayesian_t3$pcor_adj,
    labels = nodes,
    layout = "circle",
    theme = "colorblind",
    title = "Time Point 3",
    groups = item_groups_names,
    nodeNames = long_names,
    legend.cex=.55
)

# Reset the plot area.
layout(1)

# #endregion


# #region Visually comparing all graphs.

# Split the plot area.
layout(matrix(1:6, nrow = 2, byrow = TRUE))

# Plot the regularized network for the first time point.
plot(plot_network_regularization_t1)

# Plot the regularized network for the second time point.
plot(plot_network_regularization_t2)

# Plot the regularized network for the third time point.
plot(plot_network_regularization_t3)

# Plot the Bayesian graph for the first time point.
plot(plot_network_bayesian_t1)

# Plot the Bayesian graph for the second time point.
plot(plot_network_bayesian_t2)

# Plot the Bayesian graph for the third time point.
plot(plot_network_bayesian_t3)

# Reset the plot area.
layout(1)

# #endregion

## region descriptives of network structures estimated with BGGM or EBIC g-LASSO

#Most central nodes 
centralityPlot(graph_regularization_t1, scale = "raw0", include = "Strength") 
centralityPlot(graph_regularization_t2, scale = "raw0", include = "Strength")
centralityPlot(graph_regularization_t3, scale = "raw0", include = "Strength")
centralityPlot(graph_bayesian_t1$pcor_adj, scale = "raw0", include = "Strength")
centralityPlot(graph_bayesian_t2$pcor_adj, scale = "raw0", include = "Strength")
centralityPlot(graph_bayesian_t3$pcor_adj, scale = "raw0", include = "Strength")

#Global strength
global_strength_r1 <- sum(graph_regularization_t1)/2
global_strength_r1
global_strength_r2 <- sum(graph_regularization_t2)/2
global_strength_r2
global_strength_r3 <- sum(graph_regularization_t3)/2
global_strength_r3

global_strength_b1 <- sum(graph_bayesian_t1$pcor_adj)/2
global_strength_b1
global_strength_b2 <- sum(graph_bayesian_t2$pcor_adj)/2
global_strength_b2
global_strength_b3 <- sum(graph_bayesian_t3$pcor_adj)/2
global_strength_b3

#Edge weights similarity
weights_regularization_t1 <- graph_regularization_t1[lower.tri(graph_regularization_t1)]
weights_regularization_t2 <- graph_regularization_t2[lower.tri(graph_regularization_t2)]
weights_regularization_t3 <- graph_regularization_t3[lower.tri(graph_regularization_t3)]
cor(weights_regularization_t1, weights_regularization_t2)
cor(weights_regularization_t2,weights_regularization_t3)

weights_bayesian_t1 <- graph_bayesian_t1$pcor_adj[lower.tri(graph_bayesian_t1$pcor_adj)]
weights_bayesian_t2 <- graph_bayesian_t2$pcor_adj[lower.tri(graph_bayesian_t2$pcor_adj)]
weights_bayesian_t3 <- graph_bayesian_t3$pcor_adj[lower.tri(graph_bayesian_t3$pcor_adj)]
cor(weights_bayesian_t1,weights_bayesian_t2)
cor(weights_bayesian_t2,weights_bayesian_t3)

#Mean edge weight
mean(weights_regularization_t1)
mean(weights_regularization_t2)
mean(weights_regularization_t3)
mean(weights_bayesian_t1)
mean(weights_bayesian_t2)
mean(weights_bayesian_t3)
#Density
mean(weights_regularization_t1!=0)
mean(weights_regularization_t2!=0)
mean(weights_regularization_t3!=0)
mean(weights_bayesian_t1!=0)
mean(weights_bayesian_t2!=0)
mean(weights_bayesian_t3!=0)

#Number of edges estimated to be 0
sum(weights_regularization_t1==0) 
sum(weights_regularization_t2==0) 
sum(weights_regularization_t3==0) 
sum(weights_bayesian_t1==0)
sum(weights_bayesian_t2==0)
sum(weights_bayesian_t3==0)

#Centrality indicators similarity
strength_centrality_regularization_t1 <-centrality(graph_regularization_t1)$InDegree
strength_centrality_regularization_t2<-centrality(graph_regularization_t2)$InDegree
strength_centrality_regularization_t3<-centrality(graph_regularization_t3)$InDegree
cor(strength_centrality_regularization_t1, strength_centrality_regularization_t2)
cor(strength_centrality_regularization_t2, strength_centrality_regularization_t3)

strength_centrality_bayesian_t1 <-centrality(graph_bayesian_t1$pcor_adj)$InDegree
strength_centrality_bayesian_t2 <-centrality(graph_bayesian_t2$pcor_adj)$InDegree
strength_centrality_bayesian_t3 <-centrality(graph_bayesian_t3$pcor_adj)$InDegree
cor(strength_centrality_bayesian_t1,strength_centrality_bayesian_t2)
cor(strength_centrality_bayesian_t2,strength_centrality_regularization_t3)

##end region

# #region Edge weights uncertainty for the regularized approach.

# Note that `bootnet` errors on my computer if more than 1 core is used.

# Edge weights confidence intervals for the first time point.
ci_network_regularization_t1 <- bootnet(
    fit_network_regularization_t1,
    nBoots = 1000,
    nCores = 1
)

# Edge weights confidence intervals for the second time point.
ci_network_regularization_t2 <- bootnet(
    fit_network_regularization_t2,
    nBoots = 1000,
    nCores = 1
)

# Edge weights confidence intervals for the third time point.
ci_network_regularization_t3 <- bootnet(
    fit_network_regularization_t3,
    nBoots = 1000,
    nCores = 1
)

# Summary of edge weights confidence intervals for the first time point.
summary(ci_network_regularization_t1, statistics = "edge")

# Summary of edge weights confidence intervals for the second time point.
summary(ci_network_regularization_t2, statistics = "edge")

# Summary of edge weights confidence intervals for the third time point.
summary(ci_network_regularization_t3, statistics = "edge")

# Plot edge weight CI for the first time point
 plot(
 ci_network_regularization_t1,
 statistics = "edge",
 plot = "interval",
 order = "id"
 )
 # Plot edge weight CI for the second time point
 plot(
 ci_network_regularization_t2,
 statistics = "edge",
 plot = "interval",
 order = "id"
 )
 # Plot edge weight CI for the third time point
 plot(
 ci_network_regularization_t3,
 statistics = "edge",
 plot = "interval",
 order = "id"
 )

# #endregion


# #region Edge weights uncertainty for the Bayesian approach.

# Edge weights credible intervals for the first time point.
summary(fit_network_bayesian_t1)

# Edge weights credible intervals for the second time point.
summary(fit_network_bayesian_t2)

# Edge weights credible intervals for the third time point.
summary(fit_network_bayesian_t3)

# #endregion

# # region Centrality Stability

#       Regularization: 

# Run the case-dropping subset bootstrap to assess stability for the first time point
centrality_stability_t1 <- bootnet(
    fit_network_regularization_t1,
    boots=1000,
    nCores=1,
    type="case",
    useCommunities = item_groups_names,
    statistics=c("strength"),
    verbose= TRUE
 )
 #Correlation Stability Coefficient T1
 corStability(centrality_stability_t1)


 # Run the case-dropping subset bootstrap to assess stability for the second time point
centrality_stability_t2<-bootnet(
    fit_network_regularization_t2,
    boots=1000,
    nCores=1,
    type="case",
   # statistics=c("strength", "bridgeStrength"),
   statistics=c("strength"),
    verbose= TRUE
 )
 #Correlation Stability Coefficient T2
 corStability(centrality_stability_t2)


 # Run the case-dropping subset bootstrap to assess stability for the third time point
centrality_stability_t3<-bootnet(
    fit_network_regularization_t3,
    boots=1000,
    nCores=1,
    type="case",
    statistics=c("strength"),
    verbose= TRUE
 ) 
 #Correlation Stability Coefficient T3
 corStability(centrality_stability_t3) 


## end region


# #region Compare edge weights across time points for the regularized approach.

# Compare edge weights for the first and second time points.
comparison_network_regularization_t1_t2 <- NCT(
    fit_network_regularization_t1,
    fit_network_regularization_t2,
    test.edges = TRUE,
    it = 1000
)
#Global strenghts invariance time1 ~ time2
comparison_network_regularization_t1_t2$glstrinv.pval


# Compare edge weights for the second and third time points.
comparison_network_regularization_t2_t3 <- NCT(
    fit_network_regularization_t2,
    fit_network_regularization_t3,
    test.edges = TRUE,
    it = 1000
)
#Global strenghts invariance time2 ~ time3
comparison_network_regularization_t2_t3$glstrinv.pval


# Print significant edges differences for the first and second time points.
comparison_network_regularization_t1_t2$einv.pvals %>%
    .[.$`p-value` < 0.05, ]

# Print significant edges differences for the second and third time points.
comparison_network_regularization_t2_t3$einv.pvals %>%
    .[.$`p-value` < 0.05, ]

# #endregion


# #region Compare edge weights across time points for the Bayesian approach.

# Function to prepare data for comparison and avoid repetition.
prepare_data_comparison_bayesian <- function(data) {
    # Add one to the ordinal variables for the first time point.
    mutate(.data = data, across(!unhealthy, ~ .x + 1)) %>%

    # Convert to regular data frame.
    as.data.frame()
}

# Prepare first time point data for comparison.
data_comparison_bayesian_t1 <- data_t1 %>%
    prepare_data_comparison_bayesian()

# Prepare second time point data for comparison.
data_comparison_bayesian_t2 <- data_t2 %>%
    prepare_data_comparison_bayesian()

# Prepare third time point data for comparison.
data_comparison_bayesian_t3 <- data_t3 %>%
    prepare_data_comparison_bayesian()

# Compare edge weights for the first and second time points.
comparison_network_bayesian <- ggm_compare_estimate(
        data_comparison_bayesian_t1,
        data_comparison_bayesian_t2,
        data_comparison_bayesian_t3,
        formula = ~ unhealthy,
        type = "ordinal",
        impute = FALSE,

        # Check the documentation for why this is set as such.
        prior_sd = 0.5
    )

# Summarize the comparison.
summary(comparison_network_bayesian)

# Plot the comparison.
plot(summary(comparison_network_bayesian))

# Select the graphs of differences across time points.
comparison_graphs_bayesian <- select(comparison_network_bayesian)

# Split the plot area.
layout(matrix(1:3, nrow = 1))

# Plot the differences between the first and second time points.
plot_comparison_bayesian_t1_t2 <- qgraph(
    input = comparison_graphs_bayesian$pcor_adj[[1]],
    labels = nodes,
    layout = "circle",
    theme = "colorblind",
    title = "Time Point 1 - Time Point 2"
)



# Plot the differences between the first and third time points.
plot_comparison_bayesian_t1_t3 <- qgraph(
    input = comparison_graphs_bayesian$pcor_adj[[2]],
    labels = nodes,
    layout = "circle",
    theme = "colorblind",
    title = "Time Point 1 - Time Point 3"
)

# Plot the differences between the second and third time points.
plot_comparison_bayesian_t2_t3 <- qgraph(
    input = comparison_graphs_bayesian$pcor_adj[[3]],
    labels = nodes,
    layout = "circle",
    theme = "colorblind",
    title = "Time Point 2 - Time Point 3"
)

# Reset the plot area.
layout(1)

# #endregion

# #CP region

# Determine thresholds for the first time point.
thresholds_regularization_t1 <- cpThreshold(
    W = graph_regularization_t1,
    method = "weighted",
    k.range = c(3, 16),
    I.range = seq(0.40, 0.01, by = -0.005),
    threshold = c("largest.components.ratio", "chi")
)

# Inspect the thresholds for the first time point.
thresholds_regularization_t1

# Run the algorithm for the first time point using `k = 3` and `I = 0.115`.
cp_regularization_t1 <- cpAlgorithm(
    W = graph_regularization_t2,
    method = "weighted",
    k = 3,
    I = 0.115
)

# Determine thresholds for the second time point.
thresholds_regularization_t2 <- cpThreshold(
    W = graph_regularization_t2,
    method = "weighted",
    k.range = c(3, 16),
    I.range = seq(0.40, 0.01, by = -0.005),
    threshold = c("largest.components.ratio", "chi")
)

# Inspect the thresholds for the second time point.
thresholds_regularization_t2
View(thresholds_regularization_t2)

# Run the algorithm for the second time point using `k = 3` and `I = 0.125`.
cp_regularization_t2 <- cpAlgorithm(
    W = graph_regularization_t2,
    method = "weighted",
    k = 3,
    I = 0.125
)

# Determine thresholds for the third time point.
thresholds_regularization_t3 <- cpThreshold(
    W = graph_regularization_t3,
    method = "weighted",
    k.range = c(3, 16),
    I.range = seq(0.40, 0.01, by = -0.005),
    threshold = c("largest.components.ratio", "chi")
)

# Inspect the thresholds for the third time point.
thresholds_regularization_t3

# Run the algorithm for the third time point using `k = 3` and `I = 0.135`.
cp_regularization_t3 <- cpAlgorithm(
    W = graph_regularization_t3,
    method = "weighted",
    k = 3,
    I = 0.135
)

# Summarize the results for the first time point.
summary(cp_regularization_t1)

# Summarize the results for the second time point.
summary(cp_regularization_t2)

# Summarize the results for the third time point.
summary(cp_regularization_t3)

# #endregion


# #region Clique percolation for the Bayesian approach.

# Determine thresholds for the first time point.
thresholds_bayesian_t1 <- cpThreshold(
    W = graph_bayesian_t1$pcor_adj,
    method = "weighted",
    k.range = c(3, 16),
    I.range = seq(0.40, 0.01, by = -0.005),
    threshold = c("largest.components.ratio", "chi")
)

# Inspect the thresholds for the first time point.
View(thresholds_bayesian_t1)

# Run the algorithm for the first time point using `k = 3` and `I = 0.165`.
cp_bayesian_t1 <- cpAlgorithm(
    W = graph_bayesian_t1$pcor_adj,
    method = "weighted",
    k = 3,
    I = 0.125
)

# Determine thresholds for the second time point.
thresholds_bayesian_t2 <- cpThreshold(
    W = graph_bayesian_t2$pcor_adj,
    method = "weighted",
    k.range = c(3, 16),
    I.range = seq(0.40, 0.01, by = -0.005),
    threshold = c("largest.components.ratio", "chi")
)

# Inspect the thresholds for the second time point.
thresholds_bayesian_t2

# Run the algorithm for the second time point using `k = 3` and `I = 0.125`.
cp_bayesian_t2 <- cpAlgorithm(
    W = graph_bayesian_t2$pcor_adj,
    method = "weighted",
    k = 3,
    I = 0.125
)

# Determine thresholds for the third time point.
thresholds_bayesian_t3 <- cpThreshold(
    W = graph_bayesian_t3$pcor_adj,
    method = "weighted",
    k.range = c(3, 16),
    I.range = seq(0.40, 0.01, by = -0.005),
    threshold = c("largest.components.ratio", "chi")
)

# Inspect the thresholds for the third time point.
thresholds_bayesian_t3

# Run the algorithm for the third time point using `k = 3` and `I = 0.175`.
cp_bayesian_t3 <- cpAlgorithm(
    W = graph_bayesian_t3$pcor_adj,
    method = "weighted",
    k = 3,
    I = 0.125
)

# Summarize the results for the first time point.
summary(cp_bayesian_t1)

# Summarize the results for the second time point.
summary(cp_bayesian_t2)

# Summarize the results for the third time point.
summary(cp_bayesian_t3)

# #endregion


# #region Plot clique percolation plots for both approaches.

# Split the plot area.
layout(matrix(1:6, nrow = 2, byrow = TRUE))

# Plot the communities for the first time point for the regularized approach.
plot_cp_regularization_t1 <- cpColoredGraph(
    W = graph_regularization_t1,
    list.of.communities = cp_regularization_t1$list.of.communities.labels,
    theme = "colorblind",
    title = "Time Point 1"
)

# Plot the communities for the second time point for the regularized approach.
plot_cp_regularization_t2 <- cpColoredGraph(
    W = graph_regularization_t2,
    list.of.communities = cp_regularization_t2$list.of.communities.labels,
    theme = "colorblind",
    title = "Time Point 2"
)

# Plot the communities for the third time point for the regularized approach.
plot_cp_regularization_t3 <- cpColoredGraph(
    W = graph_regularization_t3,
    list.of.communities = cp_regularization_t3$list.of.communities.labels,
    theme = "colorblind",
    title = "Time Point 3"
)

# Plot the communities for the first time point for the Bayesian approach.
plot_cp_bayesian_t1 <- cpColoredGraph(
    W = graph_bayesian_t1$pcor_adj,
    list.of.communities = cp_bayesian_t1$list.of.communities.labels,
    labels = nodes,
    theme = "colorblind",
    title = "Time Point 1"
)

# Plot the communities for the second time point for the Bayesian approach.
plot_cp_bayesian_t2 <- cpColoredGraph(
    W = graph_bayesian_t2$pcor_adj,
    list.of.communities = cp_bayesian_t2$list.of.communities.labels,
    labels = nodes,
    theme = "colorblind",
    title = "Time Point 2"
)

# Plot the communities for the third time point for the Bayesian approach.
plot_cp_bayesian_t3 <- cpColoredGraph(
    W = graph_bayesian_t3$pcor_adj,
    list.of.communities = cp_bayesian_t3$list.of.communities.labels,
    labels = nodes,
    theme = "colorblind",
    avoid.repeated.mixed.colors = TRUE,
    title = "Time Point 3"
)

# Reset the plot area.
layout(1)

# #endregion

# # Region Bridge symptoms
# Estimate EGA t1
ega.r_t1 <- EGA(
  data = data_t1,
  plot.EGA = FALSE # No plot for CRAN checks
)
# Network loadings
net.loads(ega.r_t1)

# Estimate EGA t2
ega.r_t2 <- EGA(
  data = data_t2,
  plot.EGA = FALSE # No plot for CRAN checks
)
# Network loadings
net.loads(ega.r_t2)


# Estimate EGA t1
ega.r_t3 <- EGA(
  data = data_t3,
  plot.EGA = FALSE # No plot for CRAN checks
)
# Network loadings
net.loads(ega.r_t3)