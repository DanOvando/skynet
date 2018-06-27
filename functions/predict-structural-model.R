predict_structural_model <- function(mle_fit, data,
                                     mle_vars,
                                     mp = 0){

  cost_data <- data %>%
    select(mle_vars) %>%
    as.matrix()

  betas <- mle_fit$par[str_detect(names(mle_fit$par),'beta')]

  costs <-  cost_data %*% betas

  q <- exp(mle_fit$par[str_detect(names(mle_fit$par),'log_q')])

  d_hat = (exp(q * data$total_hours) * (costs + mp)) / data$aggregate_price;

  d_hat[d_hat <= 1e-3] <- 1e-3 / (2.0 - d_hat[d_hat <= 1e-3] /1e-3);

return(as.numeric(d_hat))
}