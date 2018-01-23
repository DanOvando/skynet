# include <TMB.hpp>


template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  if ( x >= eps ){
    return x;
  } else {
    pen += Type(0.01) * pow(x-eps,2);
    return eps/(Type(2.0)-x/eps);
  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  // data block

  DATA_MATRIX(data); // covariates

  DATA_VECTOR(log_d);

  DATA_VECTOR(effort);

  DATA_VECTOR(price);

  DATA_SCALAR(mp);

  DATA_INTEGER(n);

  Type nll; //blank storage for accumulated nll

  nll = 0;
  // parameter block

  PARAMETER_VECTOR(betas);

  PARAMETER(log_sigma);

  PARAMETER(logit_q);


  // model block

  Type fpen = 0.;

  Type sigma = exp(log_sigma);

  Type q = invlogit(logit_q);

  matrix<Type> cost = data * betas;

  vector<Type> d_hat(n);

  vector<Type> log_d_hat(n);


  for (int i = 0; i<n; i++){

  log_d_hat(i) = log(exp(cost(i)) + mp) - log(price(i) * exp(-q * effort(i)));

  d_hat(i) = exp(log_d_hat(i));

  }

  nll -= sum(dnorm(log_d, log_d_hat, sigma, true));

  // report block

  ADREPORT(log_d_hat);

  ADREPORT(d_hat);

  REPORT(log_d_hat);

  REPORT(d_hat);

  return nll;

}