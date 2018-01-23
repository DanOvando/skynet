predict_structural_model <- function(mle_fit, data,
                                     mle_vars,
                                     mp = 0){

  cost_data <- data %>%
    select(mle_vars) %>%
    as.matrix()

  betas <- mle_fit$par[str_detect(names(mle_fit$par),'beta')]

  costs <-  exp(cost_data %*% betas)

  q <- boot::inv.logit( mle_fit$par[str_detect(names(mle_fit$par),'logit_q')])

  log_d_hat <- log(costs + mp) - log(data$aggregate_price * exp(-q * data$total_hours));

return(as.numeric(log_d_hat))
}