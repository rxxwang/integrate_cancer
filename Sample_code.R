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
  x = rep(0:1, sum(n)),
  y = sample(0:1, sum(n), replace = T),
  strata = rep(1:(sum(n)/2), each = 2),
  cancer = rep(cancer_name, n)
)

df2 <- data.frame(
  x = rep(0:1, sum(n)),
  y = rbinom(sum(n), 1, prob = 0.4 + 0.2 * rep(0:1, sum(n))),
  strata = rep(1:(sum(n)/2), each = 2),
  cancer = rep(cancer_name, n)
)

result1 = run_integrations(df1, "Breast", F)
print(result1)

result2 = run_integrations(df2, "Breast", F)
print(result2)
