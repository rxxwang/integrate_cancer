# Asymmetric Integration of Various Cancer Datasets for Identifying Risk-Associated Variants and Genes

This project using an asymmetric integration methods for conditional logistic regression of SNPs and cancer. Here we only provide the codes for computing integrated p-values. Further codes are available upon request.

## Requirements
- R
- RcppArmadillo
- survival
- dplyr

## Files
- [infcts.condlogistic.R](./infcts.condlogistic.R): Codes for the integration method.
- [run_integrations.R](./run_integrations.R): Computing both non-Integrated p-values and Integrated p-values.
- [Sample_code.R](./Sample_code.R): An example of simulated genotype data was shown in this code.
