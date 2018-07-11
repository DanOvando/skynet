# include <TMB.hpp>


template<class Type>
Type posfun(Type x, Type eps)
{
  if ( x >= eps ){
    return x;
  } else {
    return eps/(Type(2.0)-x/eps);
  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  // data block

  DATA_MATRIX(data); // covariates

  DATA_MATRIX(location); // location covariates

  DATA_VECTOR(log_effort);

  DATA_VECTOR(price);

  DATA_SCALAR(mp);

  DATA_INTEGER(n);

  Type nll; //blank storage for accumulated nll

  nll = 0;
  // parameter block

  PARAMETER_VECTOR(betas);

  PARAMETER_VECTOR(cpue_betas);

  PARAMETER(log_mean_cpue);

  PARAMETER(log_sigma_cpue);

  PARAMETER(log_sigma);

  PARAMETER(log_q);

  // model block

  Type q = exp(log_q);

  Type sigma = exp(log_sigma);

  Type sigma_cpue = exp(log_sigma_cpue);

  Type mean_cpue = exp(log_mean_cpue);

  matrix<Type> cost = data * betas;

  matrix<Type> cpue = location * cpue_betas;

  vector<Type> effort_hat(n);

  vector<Type> log_effort_hat(n);

  for (int i = 0; i<n; i++){

  // cost(i) = posfun(cost(i), Type(0.00001), fpen);

  // d_hat(i) = posfun((exp(q * effort(i)) * (cost(i) + mp)) / price(i),
  //       Type(0.001));


    effort_hat(i) = (1 / q) * log((price(i) * cpue(i)) / (cost(i) + mp));

    //d_hat(i) = (exp(q * effort(i)) * (cost(i) + mp)) / price(i);

    if (effort_hat(i) <= 1e-3){

      effort_hat(i) = 1e-3 / (Type(2.0) - effort_hat(i) /1e-3);
    }

  log_effort_hat(i) = log(effort_hat(i));

  }

  nll -= sum(dnorm(log_effort, log_effort_hat, sigma, true));

  nll -= sum(dnorm(cpue_betas, mean_cpue, sigma_cpue, true));


  // report block

  REPORT(log_effort_hat);

  REPORT(effort_hat);

  REPORT(sigma);

  return nll;
}
