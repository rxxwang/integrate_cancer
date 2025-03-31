# Asymmetric Integration of Various Cancer Datasets for Identifying Risk-Associated Variants and Genes

This project using an asymmetric integration method for conditional logistic regression of SNPs and cancer. Here we only provide the codes for computing integrated p-values. Further codes are available upon request.

## Requirements
- R
- RcppArmadillo
- survival
- dplyr

## Files
- [infcts.condlogistic.R](./infcts.condlogistic.R): Codes for the integration method.
- [run_integrations.R](./run_integrations.R): Computing both non-Integrated p-values and Integrated p-values.
- [Sample_code.R](./Sample_code.R): An example of simulated genotype data was shown in this code.

## Data requirement
A data frame with the following features:
- x: Genotype (=0, 1, 2)
- y: Case-control status (=0, 1)
- strata: Group number of the case and control
- cancer: Cancers indicate different datasets. One is the local dataset and the others are external datasets 

## Function 
The key function is in [run_integrations.R](./run_integrations.R): **run_integrations(data, local, permuted)**
- data: A data frame with the structure described in the Data requirement
- local: One cancer name among the cancer variables of the input data frame which will be treated as the local dataset and others as external datasets.
- permuted: Decide whether to swap the case-control status (y) in order to create permuted/non-permuted p-values.

## Sample code
Shown in [Sample_code.R](./Sample_code.R).
