// SEM MNAR, yo and yu jointly sampling and user defined log lokeligood

data {
  int<lower=0> n;                // Number of data points
  int<lower=0> p;                // Number of predictors+1, regression 
  int<lower=0> q;                // Number logistic paras
  int<lower=0> no;
  int<lower=0> nu;
  
  matrix[n, (q-1)] xm;
  matrix[n, p] x;                // Predictor matrix
  matrix[nu,p] xu;
  matrix[no,p] xo;
  vector[no] yo;                   // Response vector
  int<lower=0, upper=1> m[n];
  
  matrix[n, n] w;
  matrix[n, n] I;
  matrix[n, n] wpluwt;
  matrix[n, n] wtw;
  
  vector[p] mu_beta;
  matrix[p, p] sigma2_beta;
  vector[q] mu_psi;
  matrix[q, q] sigma2_psi;
}


parameters {
  
  vector[p] beta;
  real gama;
  real lamda;
  vector[q] psi;
  vector[nu] yu; //latent (yu) treats as parameters, 
  //so should be defined as a parameter
}


model{
  // Prior distributions for parameters
  target+=multi_normal_lpdf(beta | mu_beta, sigma2_beta);    // Prior for regression coefficients
  target+=normal_lpdf(gama | 0,10000);                     
  target+=normal_lpdf(lamda |  0,10000);                         // Prior for rho
  target+=multi_normal_lpdf(psi | mu_psi, sigma2_psi);
  
  real sigma2;
  sigma2=exp(gama); // retransform
  real rho;
  rho=(exp(lamda)-1)/(exp(lamda)+1); // retransform
  
  matrix[n,n] M=I - rho * wpluwt + rho * rho * wtw;
  
  
  vector[n] y = append_row(yo, yu);
  
  target+=-0.5*n*log(2*pi())-0.5*n*log(sigma2)+
    0.5*2*sum(log(diagonal(cholesky_decompose(M))))-
    0.5*(transpose(y-x*beta)*M*(y-x*beta))/sigma2;
  
  
  
  target+=bernoulli_logit_lpmf(m | xm*psi[1:(q-1)]+y*psi[q]); // missing model
  
  
  
}

