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

  PARAMETER(log_q);

  // model block

  // Type fpen = 0.;

  Type q = exp(log_q);

  Type sigma = exp(log_sigma);

  matrix<Type> cost = data * betas;

  vector<Type> d_hat(n);

  vector<Type> log_d_hat(n);

  for (int i = 0; i<n; i++){

  // cost(i) = posfun(cost(i), Type(0.00001), fpen);

  // d_hat(i) = posfun((exp(q * effort(i)) * (cost(i) + mp)) / price(i),
  //       Type(0.001));

    d_hat(i) = (exp(q * effort(i)) * (cost(i) + mp)) / price(i);

    if (d_hat(i) <= 1e-3){

      d_hat(i) = 1e-3 / (Type(2.0) - d_hat(i) /1e-3);

      std::cout << "WTF" << '\n';

    }

  log_d_hat(i) = log(d_hat(i));

  }

  nll -= sum(dnorm(log_d, log_d_hat, sigma, true));


  // report block

  // ADREPORT(log_d_hat);

  // ADREPORT(d_hat);

  REPORT(log_d_hat);

  REPORT(d_hat);

  REPORT(sigma);

  return nll;
}
