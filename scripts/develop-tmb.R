library(TMB)
compile(here::here('scripts','fit_structural_skynet.cpp'))
dyn.load(dynlib(here::here("scripts","fit_structural_skynet")))


struct_data <- list(
  data = as.matrix(independent_data %>%
                     select(dist_from_port,
                            mean_vessel_length,
                            total_engine_power,
                            no_take_mpa,
                            restricted_use_mpa,
                            m_below_sea_level)),
  log_d = as.numeric(independent_data$log_density),
  effort = as.numeric(independent_data$total_hours),
  price = as.numeric(independent_data$aggregate_price),
  mp = 0,
  n = nrow(independent_data)
)

struct_params <- list(
  betas = rep(0, ncol(struct_data$data)),
  log_sigma = log(sd(struct_data$log_d)),
  logit_q = log(.001 / (1 - .001))
)

model <- MakeADFun(data=struct_data,parameters=struct_params)

fit <-
  nlminb(
    model$par,
    objective = model$fn,
    gradient = model$gr,
    control = list("trace" = 1)
  )

fit <-
  TMBhelper::Optimize(obj = model,
    control = list("trace" = 1)
  )
fit_report <- model$report()

check <- independent_data %>%
  mutate(log_density_hat = fit_report$log_d_hat)

check %>%
  ggplot(aes(pred, log_density_hat, color = total_hours > 0)) +
  geom_point()

