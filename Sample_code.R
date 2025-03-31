source("run_integrations.R")
source("infcts.condlogistic.R")

library(survival)
library(dplyr)

set.seed(123)  # For reproducibility

cancer_name = c("Bladder", "Brain", "Breast", "Esoph", "HN",
                "Kidney", "Leukemia", "Liver", "Lung", "Ovarian",
                "Panc", "Prostate", "Sarcoma", "Stomach")

# Generate sample data
n <- 2*sample(1:1000, 14)  # Number of observations for each cancer
df1 <- data.frame(
  x = sample(0:1, sum(n), replace = T),  # Genotype
  y = rep(0:1, sum(n)),  # Case-control status
  strata = rep(1:(sum(n)/2), each = 2),  # paired numbers of case and control
  cancer = rep(cancer_name, n)  # cancer name (local or external dataset)
)

df2 <- data.frame(
  x = rbinom(sum(n), 1, prob = 0.4 + 0.2 * rep(0:1, sum(n))),  # Generate genotype that is related with case-control status
  y = rep(0:1, sum(n)),  # Case-control status
  strata = rep(1:(sum(n)/2), each = 2),  # paired numbers of case and control
  cancer = rep(cancer_name, n)  # cancer name (local or external dataset)
)

result1 = run_integrations(df1, "Breast", F)
print(result1)

result2 = run_integrations(df2, "Breast", F)
print(result2)
