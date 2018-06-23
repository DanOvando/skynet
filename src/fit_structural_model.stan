data {
int <lower = 0> n; // number of observations

int<lower = 0> n_betas; // number of cost covariates

matrix[n, n_betas] cost_data;

vector[n] log_d;

vector[n] effort;

vector[n] price;

real mp;

real max_q;


int <lower = 0> test_n; // number of observations

matrix[test_n, n_betas] test_cost_data;

vector[test_n] test_effort;

vector[test_n] test_price;

}

parameters{

  vector<lower = 0>[n_betas] betas; //cost n_betas

  real<lower = 0> sigma; //standard deviation

  real<lower = 1e-8,upper = max_q> q;

}

transformed parameters{

vector[n] d_hat;

vector[n] log_d_hat;

vector[n] cost;

cost = cost_data * betas;

d_hat = ((exp(q * effort) .* (cost + mp)) ./ (price));

print(max(d_hat))

log_d_hat = log(d_hat);

}

model{



log_d ~ normal(log_d_hat, sigma);

betas ~ normal(0,10);

sigma ~ cauchy(0,2.5);

}

generated quantities{

vector[n] pp_d_hat;

vector[test_n] test_d_hat;

vector[test_n] mean_test_d_hat;

vector[test_n] test_cost;

test_cost = test_cost_data * betas;

mean_test_d_hat = ((exp(q * test_effort) .* (test_cost + mp)) ./ (test_price));

for (i in 1:test_n){

    test_d_hat[i] = exp(normal_rng(log(mean_test_d_hat[i]), sigma)); // pp for testing
}


for (i in 1:n){

    pp_d_hat[i] = exp(normal_rng(log_d_hat[i], sigma)); // pp for training
}


}
