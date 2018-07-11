data {
int <lower = 0> n; // number of observations

int<lower = 0> n_betas; // number of cost covariates

int<lower = 0> n_cpue_betas; // number of cost covariates

matrix[n, n_betas] cost_data;

matrix[n, n_cpue_betas] location_data;

vector[n] log_effort;

vector[n] price;

real mp;

real max_q;

}

parameters{

  vector<lower = 0> [n_betas] betas; //cost n_betas

  vector[n_cpue_betas] cpue_betas;

  real<lower = 0> mean_cpue;

  real<lower = 0> sigma; //standard deviation

  real<lower = 0> sigma_cpue;

  real<lower = .01*max_q,upper = max_q> q;

}

transformed parameters{

vector[n] effort_hat;

vector[n] log_effort_hat;

vector[n] cost;

vector[n] cpue;

vector[n] knot_cpue;

knot_cpue = mean_cpue + sigma_cpue *  cpue_betas;

cost = cost_data * betas;

cpue = location_data * knot_cpue;
print("hello")
for (i in 1:n){

  // effort_hat[i] = (1 / q) * log((price[i] * q * cpue[i]) / (cost[i] + mp));

  effort_hat[i] = (price[i] * q * cpue[i] - mp) / (2 * cost[i]);


  // (1 / q) * log((price[i] * q * cpue[i]) / (cost[i] + mp));


  if (effort_hat[i] <= 1e-3){

    effort_hat[i] =  1e-3/(2.0- effort_hat[i] / 1e-3);

  }

}

log_effort_hat = log(effort_hat);

}

model{

log_effort ~ normal(log_effort_hat, sigma);

betas ~ normal(0,10);

cpue_betas ~ normal(0, 1);

sigma ~ cauchy(0,2.5);

sigma_cpue ~ cauchy(0,2.5);

mean_cpue ~ normal(40, 10);

}
