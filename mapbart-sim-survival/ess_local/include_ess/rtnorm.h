
double rtnorm(double mean, double tau, double sd, rn& gen)
{
  double x, z, lambda;

  /* Christian Robert's way */
  //assert(mean < tau); //assert unnecessary: Rodney's way
  tau = (tau - mean)/sd;

  /* originally, the function did not draw this case */
  /* added case: Rodney's way */
  if(tau<=0.) {
    /* draw until we get one in the right half-plane */
    do { z=gen.normal(); } while (z < tau);
  }
  else {
    /* optimal exponential rate parameter */
    lambda = 0.5*(tau + sqrt(tau*tau + 4.0));

    /* do the rejection sampling */
    do {
      z = gen.exp()/lambda + tau;
      //z = lambda*gen.exp() + tau;
    } while (gen.uniform() > exp(-0.5*pow(z - lambda, 2.)));
  }

  /* put x back on the right scale */
  x = z*sd + mean;

  //assert(x > 0); //assert unnecessary: Rodney's way
  return(x);

}
