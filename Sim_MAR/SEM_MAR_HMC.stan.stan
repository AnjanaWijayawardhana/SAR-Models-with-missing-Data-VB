// SEM MAR, yo and yu jointly sampling

data {
  int<lower=0> n;                // Number of data points
  int<lower=0> p;                // Number of predictors+1
  int<lower=0> no;
  int<lower=0> nu;
  
  matrix[n, p] x;                // Predictor matrix
  matrix[nu,p] xu;
  matrix[no,p] xo;
  vector[no] yo;                   // Response vector
  matrix[n, n] w;
  matrix[n, n] I;
  matrix[n, n] wpluwt;
  matrix[n, n] wtw;
  matrix[n, n] Min;
  
  vector[p] mu_beta;
  matrix[p, p] sigma2_beta;
}


parameters {
  
  vector[p] beta;
  real gama;
  real lamda;
  vector[nu] yu;

}



model{
  // Prior distributions for parameters
  target+=multi_normal_lpdf(beta | mu_beta, sigma2_beta);    
  target+=normal_lpdf(gama | 0,10000);                     
  target+=normal_lpdf(lamda |  0,10000);                         // Prior for rho
  
  matrix[n,n] M;
  vector[n] mu;
  
  real sigma2;
  sigma2=exp(gama); // retransform
  real rho;
  rho=(exp(lamda)-1)/(exp(lamda)+1); // retransform
    
  M=(I - rho * wpluwt + rho * rho * wtw);
  mu=x*beta;
 
  
  vector[n] y = append_row(yo, yu);
  vector[n] r=y-mu;
  // target+=-0.5*log(sigma2)+0.5*log_determinant(M)-0.5*(r%*%M%*%r)/sigma2;
  
  target += -0.5 *n* log(sigma2) + 0.5 * log_determinant(M) - 0.5 * dot_product(r, M * r) / sigma2;
}

