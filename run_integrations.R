run_integrations = function(data, local, permuted) {
  cancer = sort(unique(data$cancer))
  local_num = which(cancer == local)
  
  cancer_name = c("Bladder", "Brain", "Breast", "Esoph", "HN",
                  "Kidney", "Leukemia", "Liver", "Lung", "Ovarian",
                  "Panc", "Prostate", "Sarcoma", "Stomach")
  
  # Randomly swap the case-control label
  if(permuted){
    data = data %>%
      group_by(strata) %>%
      mutate(
        swap = rbinom(1, 1, 0.5),    # For each MatchID, decide swap (0 = no swap, 1 = swap)
        y = ifelse(swap == 1, 1 - y, y)  # If swap == 1, flip 0 to 1 and 1 to 0
      ) %>%
      dplyr::select(-swap)
  }
  data = data %>% arrange(strata, desc(y))
  if (!(local %in% cancer_name)) {
    stop("Not a valid Cancer Type - make sure you are using a capital letter to start")
  }
  data_vector = rep(0, 6)
  
  #re-arrange the datasets
  datasets <- split(data[, c("x", "y", "strata")], data$cancer)
  names(datasets) <- cancer_name
  #datasets <- datasets[order(names(datasets))]
  datasets = c(datasets[local_num], datasets[-local_num])
  datasets <- lapply(datasets, function(sublist) {
    sublist$x <- matrix(sublist$x, ncol = 1)  # Converting to a column matrix
    return(list(
      x = sublist$x,
      y = sublist$y,
      strata = sublist$strata
    ))
  })
    
  result <- tryCatch({
    local_fit = clogit(datasets[[local]]$y ~ datasets[[local]]$x + strata(datasets[[local]]$strata))
    data_vector[1] = summary(local_fit)$coefficients[1]
    data_vector[2] = summary(local_fit)$coefficients[3]
    data_vector[3] = summary(local_fit)$coefficients[5]
    
    if (data_vector[3] <= 0.05) {
      integrated_fit = joint.fitting(datasets = datasets, method = "testingopt")
      weights = integrated_fit$weights
      data_vector[4] = integrated_fit$par

      #Using the weights to fit a conditional logistic model to find standard errors
      #clogit does not take 0 weights, so we only combine data that had non-zero weight
      nonzero_datasets = c(1, which(integrated_fit$weights != 0) + 1)
      
      x_all = do.call(rbind, lapply(datasets[nonzero_datasets], function(df)
        df$x))
      y_all = unname(unlist(lapply(datasets[nonzero_datasets], function(df)
        unlist(df$y))))
      strata_all = unname(unlist(lapply(datasets[nonzero_datasets], function(df)
        unlist(df$strata))))
      weights_all = rep(c(1, integrated_fit$weights[which(integrated_fit$weights != 0)]),
                        times = sapply(datasets[nonzero_datasets], function(df)
                          length(df$x)))
      
      
      x_all = as.matrix(x_all, ncol = 1)
      
      model_integrated = clogit(y_all ~ x_all + strata(strata_all),
                                method = "approximate",
                                weights = weights_all)
      
      
      if (sum(integrated_fit$weights < 1 &
              integrated_fit$weights > 0) > 0) {
        data_vector[5] = summary(model_integrated)$coefficients[4]
        data_vector[6] = summary(model_integrated)$coefficients[6]
      } else{
        data_vector[5] = summary(model_integrated)$coefficients[3]
        data_vector[6] = summary(model_integrated)$coefficients[5]
      }
    }
    T
  }, error = function(err) {
    message(paste(local, " data can not be calculated"))
    return(F)
  })
  if (!result) {
    T
  }
  return(data_vector)
}