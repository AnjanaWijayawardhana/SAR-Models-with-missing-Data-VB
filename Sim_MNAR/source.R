# Load required libraries

library(igraph)
library(MASS)
library(spdep)
library(tictoc)
library(rstan)
library(Matrix)
library(mvtnorm)
library(coda)
library(ggplot2)
library(mvnfast) # fast multivariate normal sampling
library(patchwork)
library("vctrs")
library(dplyr)
library(tidyr)
library(reshape2)


######################### Supporting Functions

#Data simulation functions

my_splitter_all_xm=function(x,xm,y,w,psi){
  n=nrow(y)
  x<-matrix(x,nrow = n)
  xm<-matrix(xm,nrow = n)
  lead_1=rep(1,n)
  z=cbind(lead_1,xm,y)%*%psi 
  pr= 1/(1+exp(-z)) 
  #m=rbinom(nrow(y),1,pr) 
  
  u=runif(nrow(y))
  
  m=as.numeric(pr>u)
  
  y_us=rep(0,n)
  
  for(i in 1:n){
    if(m[i]==1){
      y_us[i]=NA
    }else{
      y_us[i]=y[i]
    }
  }
  
  ####
  
  nu<-sum(is.na(y_us))
  no<-n-nu
  
  
  missingindexes<-which(is.na(y_us))
  xu<-x[missingindexes,]
  xo<-x[setdiff(1:n,missingindexes),]
  yo<-y[setdiff(1:n,missingindexes)]
  yu<-y[missingindexes]
  xmu<-xm[missingindexes,]
  xmo<-xm[setdiff(1:n,missingindexes),]
  u<-missingindexes
  o<-setdiff(1:n,u)
  
  w<-w[c(o,u),c(o,u)]
  mo<-m[o]
  mu_<-m[u] #To overcome the confusion with mean of prior, which is also mu
  m<-c(mo,mu_)
  
  x=rbind(matrix(xo,nrow=no),matrix(xu,nrow=nu))
  xm=rbind(matrix(xmo,nrow=no),matrix(xmu,nrow=nu))
  return(list("m"=m,"X"=x,"xm"=xm,"yo"=yo,"yu"=yu,"W"=w))
  
}  #spliting and re arranging, output, w=(woo,wou,
#wuo,wuu), x=(xo,xu),yo,yu, m=(mo,mu)

simulateSEM<-function(para,p,weightmat){
  
  
  weightmat<-as(weightmat,"CsparseMatrix")
  sigma2=para[p+3]
  rho=para[p+2]
  n=ncol(weightmat)
  m=0
  std=1
  X<-matrix(rep(c(1),n))
  
  for(i in 1:p){
    X<-cbind(X,rnorm(n,m,std))
  }
  
  b<-para[0:p+1]
  errors_part<-rnorm(n,sd=sqrt(sigma2))
  
  Y<-X%*%b+solve(Diagonal(n)-rho*weightmat,errors_part)
  
  return(list("Independent"=X[,2:(p+1)],"Dependent"=Y))
  
} # Generate SSEM, efficient generator


make.SATN.data<-function(x,xm,w,m,yo){
  w<-as.matrix(w)
  n<-ncol(w)
  nu<-sum(m)
  no<-n-nu
  I<-diag(n)
  
  
  
  x_<-cbind(1,x)
  xm<-cbind(1,xm)
  q<-ncol(xm)+1
  xo=x_[1:no,]
  xu=x_[-(1:no),]
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  p=ncol(x_)
  
  ############## priors
  sigma2_beta=1000*diag(p)
  mu_beta=rep(0,p)
  sigma2_psi=1000*diag(q)
  mu_psi=rep(0,q)
  
  stan.data.list<-list(n=n,no=no,nu=nu,xm=xm,xo=xo,xu=xu,yo=yo,x=x_,
                       p=p,q=q,w=w,I=I,wpluwt=wpluswt,m=m,
                       wtw=wtw, sigma2_beta=sigma2_beta,mu_beta=mu_beta,
                       sigma2_psi=sigma2_psi,mu_psi=mu_psi)
  
  return(stan.data.list=stan.data.list)
  
  
} # STAN data list generating functuion


# Functions to simulate from the posterior

simulate_VB_yu.theta<-function(mu,B,x,d){
  
  
  p<-ncol(B)
  m<-length(d) # missing vals, betas(col.x+1), rho+sigma2y+  psis
  
  epsilon<-rnorm(m)
  z<-rnorm(p)
  
  theta<-mu+B%*%z+as.vector(d)*epsilon
  
  #rho
  lamda<-theta[nu+ncol(x)+1+1+1]
  rho<-(exp(lamda) -
          1)/(exp(lamda) +1)
  
  #sigma2y
  sigma2y<-exp(theta[nu+ncol(x)+1+1])
  
  theta[nu+ncol(x)+1+1+1]<-rho
  theta[nu+ncol(x)+1+1]<-sigma2y
  return(as.vector(theta))
  
  #################
  yu<-theta_u[1:nu]
  beta<-theta_u[(nu+1):(nu+p)]
  gama<-theta_u[(nu+p+1)]
  lamda<-theta_u[(nu+p+2)]
  psi<-theta_u[(nu+p+2+1):length(theta_u)]
  rho<-(exp(lamda)-1)/(exp(lamda)+1)
  
  
  
} # Simulate draws from JVB

simulate_VB_theta<-function(mu,B,x,d){
  
  
  p<-ncol(B)
  m<-length(d)
  
  epsilon<-rnorm(m)
  z<-rnorm(p)
  
  theta<-mu+B%*%z+as.vector(d)*epsilon
  
  #rho
  lamda<-theta[ncol(x)+1+1+1]
  rho<-(exp(lamda) -
          1)/(exp(lamda) +1)
  
  #sigma2y
  sigma2y<-exp(theta[ncol(x)+1+1])
  
  theta[ncol(x)+1+1+1]<-rho
  theta[ncol(x)+1+1]<-sigma2y
  return(as.vector(theta))
  
  
  
  
} # Simulate draws from HVB  



# Functions to find starting values

estimate_regression<-function(x,yo,no){
  # x<-cbind(1,x)
  fit.regression<-lm(yo~x[1:no,])
  return(fit.regression)
}



sim_yugyotheta_<-function(yo,listw_x_m,    
                          yu.start,rc,paras,
                          N1,k,bsize,reminder){
  nu=length(yu.start)
  n<-ncol(listw_x_m$w.list[[1]])
  I<-Diagonal(n)
  no<-length(yo)
  nu<-n-no
  ####
  
  beta<-paras[(1:rc)]
  gamma.lamda<-paras[(rc+1):(rc+2)]
  rho<-gamma.lamda[2]
  sigma2y<-(gamma.lamda[1])
  psi<-paras[-(1:(rc+2))]
  ####
  
  cal.con.dis<-function(x,yo,yui_,w,I,beta,rho,sigma2y,bsize,no,nu,reminder,j,k){ #j is the iteration number
    
    M<-t(I-rho*w)%*%(I-rho*w)
    # Sigma<-sigma2y*solve(t(I-rho*w)%*%(I-rho*w))
    if(reminder==0){
      # Sigma_ui.ui_o<-Sigma[1:bsize,(bsize+1):n] #between cov(yui and yui_,yo)
      # Sigma_ui_o.ui_o<-Sigma[(bsize+1):n,(bsize+1):n] #cov(yui_,yo and yui_,yo)
      # Sigma_uiui<-Sigma[1:bsize,1:bsize]
      # 
      # 
      # mu_uigui_o<-x[1:bsize,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
      #                                                      c(yui_,yo)-(x[(bsize+1):n,])%*%beta)
      # sigma_uigui_o<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
      
      
      # sigma_uj.v<-Sigma[1:bsize,(bsize+1):n]
      # sigma_v.v<-Sigma[(bsize+1):n,(bsize+1):n] 
      # sigma_uj.uj<-Sigma[1:bsize,1:bsize] 
      # mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
      # v<-c(yui_,yo)
      # 
      # mu_ujguv<-x[1:bsize,]%*%beta+sigma_uj.v%*%solve(sigma_v.v,(v-mu_v))
      # sigma_ujgv<-sigma_uj.uj-sigma_uj.v%*%solve(sigma_v.v,t(sigma_uj.v))
      #######################################
      M_v.v=M[(bsize+1):n,(bsize+1):n]
      M_uj.uj=M[1:bsize,1:bsize]
      M_uj.v=M[1:bsize,(bsize+1):n]
      M_v.uj=t(M_uj.v)
      v<-c(yui_,yo)
      mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
      
      In<-Diagonal(ncol(M_uj.uj))
      cholMuj.uj<-Cholesky(M_uj.uj)
      mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
      sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
      
      
      # j=3
    }else{
      if(j<=k){
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(yui_,yo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        # dim(t(Sigma_ui.ui_o))
        
      }else{ # for the falf block
        # as.matrix(Sigma)[1:reminder,(reminder+1):n,drop = FALSE]
        # Sigma_ui.ui_o<-Sigma[1:reminder,(reminder+1):n,drop = FALSE] #between cov(yui and yui_,yo)
        # Sigma_ui_o.ui_o<-Sigma[(reminder+1):n,(reminder+1):n] #cov(yui_,yo and yui_,yo)
        # Sigma_uiui<-Sigma[1:reminder,1:reminder] 
        # 
        # mu_ujguv<-x[1:reminder,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
        #                                                       c(yui_,yo)-(x[-(1:reminder),])%*%beta)
        # sigma_ujgv<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
        
        ################################
        
        M_v.v=M[(reminder+1):n,(reminder+1):n]
        M_uj.uj=M[1:reminder,1:reminder]
        M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
        M_v.uj=t(M_uj.v)
        v<-c(yui_,yo)
        mu_v<-c((x[-(1:reminder),])%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        
        
      }
      
    }
    return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
  }
  log.target<-function(m,yo,yu,x,psi){
    
    y=c(yu,yo)
    z=cbind(x,y)%*%psi
    pxy<-1/(1+exp(-z))
    
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    return(loglike1)
  }
  
  
  yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
  yu.chain[1,]<-yu.start
  current_state<-yu.start
  
  # tic()
  accepts<-0
  for(i in 2:N1){#i=1 N=5
    
    
    # j=1
    current<-0
    for (j in 1:(k+ifelse(reminder>0,1,0))) {
      # accepts<-0
      # yubj<-current_state[(current+1):(current+bsize)]
      if(j<=k){ # for full blocks j=1
        yubj_<-current_state[-((current+1):(current+bsize))]
        
        co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                            rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
        
        con.mean<-co.dis$mu_ujguv
        con.var<-co.dis$sigma_ujgv
        
        yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
        proposal_state<-current_state
        
        proposal_state[(current+1):(current+bsize)]<-yubj.star
        
        # current<-5
        ration=min(0,log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j]],psi=psi)-
                     log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j]],psi=psi))
        
        
        
        if ((!is.na(ration))){
          if((log(runif(1)) <= ration))
          {
            current_state <- proposal_state
            accepts<-accepts+1
          }
          
        }
        
        current<-current+bsize
      }else{ # for the half block j=2 current=5
        yubj_<-current_state[-((current+1):(current+reminder))]
        # current=10
        # print("in half block")
        co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                            rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
        # print("out half block")
        con.mean<-co.dis$mu_ujguv
        con.var<-co.dis$sigma_ujgv
        
        yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
        proposal_state<-current_state
        
        proposal_state[(current+1):(current+reminder)]<-yubj.star
        
        # current<-5
        ration=min(0,log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j]],psi=psi)-
                     log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j]],psi=psi))
        
        
        
        if ((!is.na(ration))){
          if((log(runif(1)) <= ration))
          {
            current_state <- proposal_state
            accepts<-accepts+1
          }
          
        }
        
      }
    }
    
    
    
    
    yu.chain[i,] <- current_state
  }
  # toc()
  
  
  
  return(list(yugyom.chain=yu.chain,accepts=accepts,
              co.dis=co.dis))
  
}

# Function to rearrange matrices W, X, and y to update y_u block-wise.

make.w_x_<-function(w,x,xm,m,k,reminder){
  
  if(reminder==0){
    
    w.list<-list()
    x.list<-list()
    xm.list<-list()
    m.list<-list()
    
    all.index<-1:n
    o<-1:no
    u<-setdiff(all.index,o)
    current<-0
    for(i in 1:k){ #i=1
      ui<-(no+1+current):(no+current+bsize)
      ui_<-setdiff(u,ui)
      
      w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
      # w[1:15,1:15]
      xui<-matrix(x[ui,],nrow=bsize,)
      xui_<-as.matrix(x[ui_,],nrow=nu-bsize)
      xo<-as.matrix(x[o,])
      x.list[[i]]<-rbind(xui,rbind(xui_,xo))
      
      xmui<-matrix(xm[ui,],nrow=bsize,)
      xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
      xmo<-as.matrix(xm[o,])
      xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
      
      m.list[[i]]<-rev(m)
      # m<-rev(m)
      current<-current+bsize
      # print("no half block")
    }}else{
      
      w.list<-list()
      x.list<-list()
      xm.list<-list()
      m.list<-list()
      
      all.index<-1:n
      o<-1:no
      u<-setdiff(all.index,o)
      current<-0
      # i=2 current=5
      for(i in 1:(k+1)){
        if(i<=k){
          ui<-(no+1+current):(no+current+bsize)
          ui_<-setdiff(u,ui)
          
          w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
          # w[1:15,1:15]
          xui<-matrix(x[ui,],nrow=bsize,)
          xui_<-matrix(x[ui_,],nrow=nu-bsize)
          xo<-x[o,]
          x.list[[i]]<-rbind(xui,rbind(xui_,xo))
          
          xmui<-matrix(xm[ui,],nrow=bsize,)
          xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
          xmo<-as.matrix(xm[o,])
          xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
          
          
          m.list[[i]]<-rev(m)
          
          current<-current+bsize
        }else{
          ui<-(no+1+current):(no+current+reminder)
          ui_<-setdiff(u,ui)
          
          w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
          # w[1:15,1:15]
          xui<-matrix(x[ui,],nrow=reminder,)
          xui_<-as.matrix(x[ui_,],nrow=nu-reminder)
          xo<-as.matrix(x[o,])
          x.list[[i]]<-rbind(xui,rbind(xui_,xo))
          
          xmui<-matrix(xm[ui,],nrow=reminder,)
          xmui_<-as.matrix(xm[ui_,],nrow=nu-reminder)
          xmo<-as.matrix(xm[o,])
          xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
          
          
          m.list[[i]]<-rev(m)
          
        }
      }
      # print("There is a half block")
    }
  
  
  return(list(w.list=w.list,x.list=x.list,xm.list=xm.list,m.list=m.list))
} 



######################### VB algorithms


MNAR_thetayu<-function(x,xm,m,w,yo,p,   # JVB algorithm
                       start_theta,N){ 
  # p=10
  w<-as(w,"CsparseMatrix")
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=ncol(xm)+1     #number of logistic covariates. i.e. betas
  totpara<-rc+lc+2
  adapt_epsilon=10^(-6)
  v=0.95
  
  n<-length(m)
  no<-length(yo)
  nu<-n-no
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  
  xo=x[1:no,]
  xu=x[-(1:no),]
  mu_<-m[-(1:no)] #To overcome the confusion with mean of prior, which is also mu
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  
  
  grad_2<-function(theta_u,wpluswt,wtw,w,x,xm,I,yo,xo,xu,m,mu_,nu,no,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p means the dimension of X.
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)] #transform of sigma2
    lamda<-theta_u[(nu+p+2)] #transform of sigma rho
    psi<-theta_u[(nu+p+3):length(theta_u)]
    
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    omega_psi=10000
    
    
    y=c(yo,yu)
    A=I-rho*w
    # M=t(A)%*%A
    M=I-rho*wpluswt+rho*rho*wtw
    #betas
    div_h_beta<-exp(-gama)*t(y-x%*%beta)%*%M%*%x-t(beta)/omega_beta
    
    #gamma
    div_h_gamma<--0.5*n+(0.5*exp(-gama))*t(y-x%*%beta)%*%M%*%(y-x%*%beta)-gama/omega_gamma
    
    #Lamda
    drho_by_dlamda<-2*exp(lamda)/(1+exp(lamda))^2
    
    dM_by_dlamda<-(-wpluswt+2*rho*wtw)*drho_by_dlamda
    # dlogdetM_by_dlamda<-sum(diag(solve(M)%*%dM_by_dlamda))
    
    # tic()
    # M.in.dMbdlamda<-solve(M,dM_by_dlamda) #1
    # print(paste("direct", sum(diag(M.in.dMbdlamda))))
    # toc()
    
    M.in.dMbdlamda<-solve(Cholesky(M),dM_by_dlamda) #2, faster
    # print(paste("Choleskey", sum(diag(M.in.dMbdlamda))))
    
    dlogdetM_by_dlamda<-sum(diag(M.in.dMbdlamda))
    
    div_h_lamda<-0.5*dlogdetM_by_dlamda-0.5*exp(-gama)*t(y-x%*%beta)%*%dM_by_dlamda%*%(y-x%*%beta)-lamda/omega_lamda
    #Psi
    z<-cbind(xm,y)
    div_h_psi<-(t(m-(exp(z%*%psi)/((1+exp(z%*%psi)))))%*%z)-t(psi)/omega_psi
    
    
    
    #yu
    
    zu<-z[(no+1):(no+nu),]
    #yu i=1
    
    psi_y<-psi[length(psi)]
    
    
    div_h_yu_1<-psi_y*(mu_-(exp(zu%*%psi))/(1+exp(zu%*%psi)))
    
    #1-This one is correct
    
    
    #2 1-This one is correct
    div_h_y<--exp(-gama)*t(y-x%*%beta)%*%M
    div_h_yu_2<-div_h_y[,(no+1):n]
    
    #3
    
    
    # #4 1-This one is correct
    div_h_yu_2<-(-exp(-gama)*t(y-x%*%beta)%*%M)[,(no+1):n]
    # dim(div_h_yu_2)
    
    div_h_yu<-div_h_yu_1+div_h_yu_2
    
    fullgrad<-c(as.vector(div_h_yu),as.vector(div_h_beta),
                as.vector(div_h_gamma),as.vector(div_h_lamda),as.vector(div_h_psi))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  
  cal_Lbound_2<-function(theta_u,nu,mu,d,B,m,xm,x,yo,wpluswt,wtw){
    
    p<-ncol(x) #length of betas
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)]
    lamda<-theta_u[(nu+p+2)]
    psi<-theta_u[(nu+p+2+1):length(theta_u)]
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    
    p_<-ncol(B)
    y<-c(yo,yu)
    n<-length(y)
    A=I-rho*w
    M=t(A)%*%A
    nu_theta<-length(theta_u)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_psi=10000
    omega_lamda=10000
    
    D<-Diagonal(n=length(d),x=as.vector(d)) #making D sparse and diagonal
    
    
    z=cbind(xm,y)%*%psi
    pxy<-1/(1+exp(-z))
    
    loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
    
    # log_det<-2*sum(log(diag(chol(M))))
    
    log_det<-2*determinant(Cholesky(M))[[1]][[1]]
    loglike2<--0.5*n*log(2*pi)-0.5*n*gama+
      0.5*log_det-0.5*exp(-gama)*t(y-x%*%beta)%*%M%*%(y-x%*%beta) # p_y efficient
    
    loglike3<-dmvnorm(x=as.vector(beta),mean = rep(0,ncol(x)),sigma=omega_beta*diag(ncol(x)),
                      log=T,checkSymmetry = F) #p_beta
    loglike4<-dnorm(x=gama,mean = 0,sd=omega_gamma) #p_gama (variance)
    
    loglike5<-dnorm(x=lamda,mean = 0,sd=omega_lamda) #p_lamda (variance)
    
    loglike6<-dmvnorm(x=as.vector(psi),mean = rep(0,(ncol(xm)+1)),sigma=omega_psi*diag(ncol(xm)+1),
                      log=T,checkSymmetry = F) #p_psi
    
    # logq<-dmvnorm(x=as.vector(theta_u),mean = as.vector(mu),sigma =B%*%t(B)+D^2,log = T,checkSymmetry = F) # q1
    
    ################################################################################ Q using Davids method
    
    r<-theta_u-mu
    B_tilda<-(1/d)*B
    r_tilde<- (1/d)*r
    
    # D_inv<-Diagonal(n=length(d),x=(1/as.vector(d)))
    # r_tilde<-D_inv%*%r #sum((1/d)*r==D_inv%*%r )
    # B_tilda<-D_inv%*%B #sum((1/d)*B ==B_tilda)
    
    
    epsilon_hat<-t(B_tilda)%*%r_tilde
    I_BtB<-forceSymmetric(Diagonal(n=p_,x=1)+t(B_tilda)%*%B_tilda) # o/w, i.e if not "forceSymmetric" logq may generate NaN values 
    R<-(Matrix::chol(I_BtB))
    # R<-t(Matrix::Cholesky(I_BtB))
    
    r_tilde2<-solve(t(R),epsilon_hat)
    rtcovr<-t(r_tilde)%*%r_tilde-t(r_tilde2)%*%r_tilde2 #Davids
    
    
    
    # rtcovr<-t(r_tilde)%*%r_tilde-t(epsilon_hat)%*%solve(I_BtB)%*%epsilon_hat #mine
    
    
    
    logdet<-sum(log(as.vector(d^2)))+sum(log(diag(R^2)))
    
    logq<--0.5*nu_theta*log(2*pi)-0.5*logdet-0.5*rtcovr #q2
    
    
    ##########################################################################
    
    # logdet<-2*sum(log(diag(chol(B%*%t(B)+D^2))))
    # logq<--0.5*nu_theta*log(2*pi)-0.5*logdet-0.5*t(r)%*%solve(B%*%t(B)+D^2,r) #q3
    
    
    
    
    return(as.numeric(loglike1)+as.numeric(loglike2)+
             as.numeric(loglike3)+as.numeric(loglike4)+
             as.numeric(loglike5)+as.numeric(loglike6)-as.numeric(logq))
    
    
  } #to calculate lower bound
  
  # 1 calculate det(M) usigng chol(), and 2 uses choleskey(). 2 is very fast
  
  
  #to calculate lower bound
  #initial values for lamdha
  #p=2
  # mu<-rep(0.1,totpara+nu)
  mu<-start_theta
  d<-rep(0.1,totpara+nu)
  B<-matrix(rep(0.1,(totpara+nu)*p),ncol=p)
  B[upper.tri(B)] <- 0
  # B<-as(B, "dgCMatrix")
  
  #initial values for adaptive learning
  E_g2_mu<-rep(0,totpara+nu)
  E_g2_d<-rep(0,totpara+nu)
  E_g2_B<-rep(0,totpara+nu) # this needs to be changed when p>1
  
  E_delta2_mu<-rep(0,totpara+nu)
  E_delta2_d<-rep(0,totpara+nu)
  E_delta2_B<-rep(0,totpara+nu) # this needs to be changed when p>1
  
  rho_mu<-rep(0,totpara+nu)
  rho_B<-rep(0,totpara+nu)
  rho_d<-rep(0,totpara+nu)
  
  # creating a vector for storing lowerbounds
  Lbound<-c()
  alllamda<-c() #To store lamda in each iteration
  tlowerbound<-c()
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) 
  #i=1
  for(i in 1:N){
    
    zeta<-(rnorm(n=(p+totpara+nu)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    z<-matrix(zeta[1:p],ncol = 1)
    
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    
    
    theta<-mu+B%*%z+d*as.vector(epsilon)
    # theta_t<-matrix(theta_t,ncol=1)
    D<-diag(d)
    
    # Calculate log gradient of theta
    
    # tic()
    gradh_theta<-grad_2(theta_u=theta,w=w,x,I=I,yo=yo,xo=xo,xu=xu,xm=xm,
                        wpluswt=wpluswt,wtw=wtw,nu=nu,no=no,n=n,mu_=mu_,m=m) # gradient of log h(theta)
    
    # toc()
    
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    
    
    #Construct UB estimates
    grad_mu<-gradh_theta+part2
    grad_B<-gradh_theta%*%t(z)+part2%*%t(z)
    gradd<-diag(gradh_theta%*%t(epsilon)+part2%*%t(epsilon))
    grad_d<-matrix(gradd,ncol = 1)
    
    #Set new learning rates
    
    # calculate new learning rates####
    
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B<-B+rho_B*grad_B
    B[upper.tri(B)]<-0
    d<-d+rho_d*grad_d
    
    # calculate lowerbound
    # tic()
    Lbound[i]<-cal_Lbound_2(theta_u=theta,nu=nu,mu=mu,d=d,B=B,m=m,xm=xm,
                            x=x,yo=yo,wpluswt=wpluswt,wtw=wtw)
    # toc()
    
    
    # tlowerbound[i]=system.time(Lbound[i]<-cal_Lbound(theta_u=theta,nu=nu,mu=mu,d=d,B=B,m=m,
    # x=x,yo=yo,wpluswt=wpluswt,wtw=wtw))[[3]]
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    
    
    alllamda[i]<-mu[nu+rc+2]
    # Lbound[i]<-0
    all_paras[i,]<-as.vector(mu)
    # print(i)
  }
  
  return(list(mu=mu,B=B,d=d,L.bound=Lbound,alllamda=
                alllamda,all_paras=all_paras))
} 



MNAR_VB_thetaAug_NoB<-function(x,xm,m,w,yo,startyu,p,   #HVB-NoB algorithm
                               start_theta,N,N1){ 
  # p=10
  w<-as(w,"CsparseMatrix")
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=ncol(xm)+1     #number of logistic covariates. i.e. betas
  totpara<-rc+lc+2
  adapt_epsilon=10^(-6)
  v=0.95
  
  n<-length(m)
  no<-length(yo)
  nu<-n-no
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  
  grad_2<-function(theta_u,wpluswt,wtw,w,x,I,yo,xo,xu,xm,m,mu_,nu,no,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p means the dimension of X.
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)] #transform of sigma2
    lamda<-theta_u[(nu+p+2)] #transform of sigma rho
    psi<-theta_u[(nu+p+3):length(theta_u)]
    
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    omega_psi=10000
    
    
    y=c(yo,yu)
    # A=I-rho*w
    # M=t(A)%*%A
    M=I-rho*wpluswt+rho*rho*wtw
    #betas
    div_h_beta<-exp(-gama)*t(y-x%*%beta)%*%M%*%x-t(beta)/omega_beta
    
    #gamma
    div_h_gamma<--0.5*n+(0.5*exp(-gama))*t(y-x%*%beta)%*%M%*%(y-x%*%beta)-gama/omega_gamma
    
    #Lamda
    drho_by_dlamda<-2*exp(lamda)/(1+exp(lamda))^2
    
    dM_by_dlamda<-(-wpluswt+2*rho*wtw)*drho_by_dlamda
    # dlogdetM_by_dlamda<-sum(diag(solve(M)%*%dM_by_dlamda))
    
    # tic()
    # M.in.dMbdlamda<-solve(M,dM_by_dlamda) #1
    # print(paste("direct", sum(diag(M.in.dMbdlamda))))
    # toc()
    
    M.in.dMbdlamda<-solve(Cholesky(M),dM_by_dlamda) #2, faster
    # print(paste("Choleskey", sum(diag(M.in.dMbdlamda))))
    
    dlogdetM_by_dlamda<-sum(diag(M.in.dMbdlamda))
    
    div_h_lamda<-0.5*dlogdetM_by_dlamda-0.5*exp(-gama)*t(y-x%*%beta)%*%dM_by_dlamda%*%(y-x%*%beta)-lamda/omega_lamda
    #Psi
    z<-cbind(xm,y)
    div_h_psi<-(t(m-(exp(z%*%psi)/((1+exp(z%*%psi)))))%*%z)-t(psi)/omega_psi
    
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_gamma),as.vector(div_h_lamda),as.vector(div_h_psi))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  sim_yugyotheta<-function(x,xu,xo,w,yo,yu,m,beta,sigma2y,rho,psi,N){
    nu=length(yu)
    samples<-matrix(rep(0,nu*N),ncol = nu)
    
    accepts<-0
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yo,yu)
      z=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-z))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    con.dis.ugo<-function(beta,rho,sigma2y,yo,xo,xu,w,I){
      A=I-rho*w
      M<-t(A)%*%A
      V<-solve(Cholesky(M),I)
      no<-nrow(xo)
      # M=I-rho*wpluswt+rho*rho*wtw
      V_oo=V[1:no,1:no]
      V_uu=V[(no+1):n,(no+1):n]
      V_ou=V[1:no,(no+1):n]
      V_uo=t(V_ou)
      
      mu_ugo=xu%*%beta+V_uo%*%solve(V_oo,(yo-xo%*%beta))
      sigma_ugo=forceSymmetric(sigma2y*(V_uu-V_uo%*%solve(V_oo,V_ou)))
      
      return(list(as.vector(mu_ugo),
                  as.matrix(sigma_ugo)))
    }
    
    current_state<-yu
    samples[1,]<-current_state
    co.dis<-con.dis.ugo(beta,rho,sigma2y,yo,xo,xu,w,I)
    # tic()
    for(i in 2:N){
      
      
      con.mean<-co.dis[[1]]
      con.var<-co.dis[[2]]
      
      # yu.star=rmvnorm(1,con.mean,con.var) #new proposal
      yu.star=rmvn(1,as.vector(con.mean),as.matrix(con.var))
      ration=min(0,log.target(m=m,yo=yo,yu=yu.star,x=x,psi=psi)-log.target(m=m,yo=yo,yu=current_state,x=x,psi=psi))
      
      
      if (log(runif(1)) <= ration) {
        current_state <- yu.star
        accepts<-accepts+1
      }
      
      
      samples[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(samples=samples,accepts=accepts))
    
  }
  
  sim_yugyotheta_<-function(x,xu,xo,w,wpluswt,wtw,yo,# efficiently calculate con mean and var
                            yu,m,beta,sigma2y,rho,psi,N){
    nu=length(yu)
    samples<-matrix(rep(0,nu*N),ncol = nu)
    
    accepts<-0
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yo,yu)
      z=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-z))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    con.dis.ugo<-function(beta,rho,sigma2y,yo,xo,xu,w,wpluswt,wtw,I){
      M=I-rho*wpluswt+rho*rho*wtw
      no<-nrow(xo)
      nu<-nrow(xu)
      # tic()
      # V<-solve(Cholesky(M),I)
      # 
      # # M=I-rho*wpluswt+rho*rho*wtw
      # V_oo=V[1:no,1:no]
      # V_uu=V[(no+1):n,(no+1):n]
      # V_ou=V[1:no,(no+1):n]
      # V_uo=t(V_ou)
      # mu_ugo=xu%*%beta+V_uo%*%solve(V_oo,(yo-xo%*%beta))
      # sigma_ugo=forceSymmetric(sigma2y*(V_uu-V_uo%*%solve(V_oo,V_ou)))
      # toc()
      
      # tic()
      M_oo=M[1:no,1:no]
      M_uu=M[(no+1):n,(no+1):n]
      M_ou=M[1:no,(no+1):n]
      M_uo=t(M_ou)
      In<-Diagonal(nu)
      cholMuu<-Cholesky(M_uu)
      mu_ugo_=xu%*%beta-solve(cholMuu,M_uo)%*%(yo-xo%*%beta)
      sigma_ugo_=forceSymmetric(sigma2y*solve(cholMuu,In))
      # toc()
      
      return(list(as.vector(mu_ugo_),
                  as.matrix(sigma_ugo_)))
    }
    
    current_state<-yu
    samples[1,]<-current_state
    # tic()
    co.dis<-con.dis.ugo(beta,rho,sigma2y,yo,xo,xu,w,wpluswt,wtw,I)
    # toc()
    # tic()
    for(i in 2:N){
      
      
      con.mean<-co.dis[[1]]
      con.var<-co.dis[[2]]
      
      # yu.star=rmvnorm(1,con.mean,con.var) #new proposal
      # tic()
      yu.star=rmvn(1,as.vector(con.mean),as.matrix(con.var))
      # toc()
      ration=min(0,log.target(m=m,yo=yo,yu=yu.star,x=x,psi=psi)-log.target(m=m,yo=yo,yu=current_state,x=x,psi=psi))
      
      
      if (log(runif(1)) <= ration) {
        current_state <- yu.star
        accepts<-accepts+1
      }
      
      
      samples[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(samples=samples,accepts=accepts))
    
  }
  
  sim_yugyotheta__<-function(x,xu,xo,xm,w,wpluswt,wtw,yo,# efficiently calculate con mean and var
                             yu,m,beta,
                             sigma2y,rho,psi,N){
    nu=length(yu)
    samples<-matrix(rep(0,nu*N),ncol = nu)
    
    accepts<-0
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yo,yu)
      z=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-z))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    con.dis.ugo<-function(beta,rho,sigma2y,yo,xo,xu,w,wpluswt,wtw,I){
      M=I-rho*wpluswt+rho*rho*wtw
      no<-nrow(xo)
      nu<-nrow(xu)
      # tic()
      # V<-solve(Cholesky(M),I)
      # 
      # # M=I-rho*wpluswt+rho*rho*wtw
      # V_oo=V[1:no,1:no]
      # V_uu=V[(no+1):n,(no+1):n]
      # V_ou=V[1:no,(no+1):n]
      # V_uo=t(V_ou)
      # mu_ugo=xu%*%beta+V_uo%*%solve(V_oo,(yo-xo%*%beta))
      # sigma_ugo=forceSymmetric(sigma2y*(V_uu-V_uo%*%solve(V_oo,V_ou)))
      # toc()
      
      # tic()
      M_oo=M[1:no,1:no]
      M_uu=M[(no+1):n,(no+1):n]
      M_ou=M[1:no,(no+1):n]
      M_uo=t(M_ou)
      In<-Diagonal(nu)
      cholMuu<-Cholesky(M_uu)
      mu_ugo_=xu%*%beta-solve(cholMuu,M_uo)%*%(yo-xo%*%beta)
      sigma_ugo_=forceSymmetric(sigma2y*solve(cholMuu,In))
      # toc()
      
      return(list(as.vector(mu_ugo_),
                  as.matrix(sigma_ugo_)))
    }
    
    
    # tic()
    co.dis<-con.dis.ugo(beta,rho,sigma2y,yo,xo,xu,w,wpluswt,wtw,I)
    # toc()
    # tic()
    con.mean<-co.dis[[1]]
    con.var<-co.dis[[2]]
    
    # tic()
    yu.star.mat=rmvn(N,as.vector(con.mean),as.matrix(con.var))
    # toc()
    yu.star.mat[1,]<-yu
    yu.chain<-matrix(rep(0,ncol(yu.star.mat)*nrow(yu.star.mat)),ncol=ncol(yu.star.mat))
    yu.chain[1,]<-yu
    current_state<-yu
    # tic()
    for(i in 2:N){#i=2
      
      # yu.star=rmvn(1,as.vector(con.mean),as.matrix(con.var))
      yu.star=yu.star.mat[i,]
      ration=min(0,log.target(m=m,yo=yo,yu=yu.star,x=xm,psi=psi)-log.target(m=m,yo=yo,yu=current_state,x=xm,psi=psi))
      
      
      if ((!is.na(ration))){
        if((log(runif(1)) <= ration))
        {
          current_state <- yu.star
          accepts<-accepts+1
        }
        
      }
      
      
      yu.chain[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(samples=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }#simultanious simulations
  
  
  #to calculate lower bound
  #initial values for lamdha
  #p=2
  #mu<-rep(0.1,totpara)
  
  mu<-c(start_theta)
  startyu<-yu
  d<-rep(0.1,totpara)
  B_vec<-matrix(rep(0.1,(totpara)*p),ncol=1)
  B<-matrix(as.vector(B_vec),ncol=p)
  B[upper.tri(B)] <- 0
  # B<-as(B, "dgCMatrix")
  
  #initial values for adaptive learning
  E_g2_mu<-rep(0,totpara)
  E_g2_d<-rep(0,totpara)
  E_g2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  E_delta2_mu<-rep(0,totpara)
  E_delta2_d<-rep(0,totpara)
  E_delta2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  rho_mu<-rep(0,totpara)
  rho_B<-rep(0,totpara*p)
  rho_d<-rep(0,totpara)
  
  # creating a vector for storing lowerbounds
  Lbound<-c()
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) #To store lamda in each iteration
  tlowerbound<-c()
  accepts<-c()
  #i=1
  for(i in 1:N){
    
    # thetas
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    z<-matrix(zeta[1:p],ncol = 1)
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    theta<-mu+B%*%z+d*as.vector(epsilon)
    D<-diag(d)
    
    # generate yus from p(yu | yo, theta)
    beta<-theta[(1:rc)]
    gamma.lamda<-theta[(rc+1):(rc+2)]
    rho<-(exp(gamma.lamda[2])-1)/(exp(gamma.lamda[2])+1)
    sigma2y<-exp(gamma.lamda[1])
    psi<-theta[-(1:(rc+2))]
    # N1=100
    
    mysamples<-sim_yugyotheta__(x=x,xu=xu,xo=xo,xm=xm,w=w,wpluswt = wpluswt,wtw =wtw ,
                                yu=startyu,yo=yo,m=m,
                                beta=beta,sigma2y=sigma2y,rho=rho,psi=psi,N=N1)
    # plot(ts((mysamples$samples)[,1]))
    accepts[i]<-mysamples$accepts
    # dim(mysamples$samples)
    # yunote<-apply((mysamples$samples)[-(1:N1/2),], 2, mean)
    yunote<-(mysamples$samples)[N1,]
    startyu<-yunote # set starting yu of the MCMC samples of the next iteration
    
    theta_u<-c(as.vector(yunote),as.vector(theta)) # combine thetas and yu
    
    # tic()
    gradg_theta<-grad_2(theta_u=theta_u,w=w,x,I=I,yo=yo,xo=xo,xu=xu,xm=xm,
                        wpluswt=wpluswt,wtw=wtw,nu=nu,no=no,n=n,mu_=mu_,m=m) # gradient of log h(theta)
    
    # toc()
    
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    
    #Construct UB estimates
    grad_mu<-gradg_theta-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradg_theta))
    grad_B<-t(t(z)%x%I_theta)%*%(gradg_theta-dlogq_by_dtheta)
    gradd<-epsilon*(gradg_theta-dlogq_by_dtheta)
    grad_d<-matrix(gradd,ncol = 1)
    
    
    #Set new learning rates
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B_vec<-B_vec+rho_B*grad_B
    B<-matrix(as.vector(B_vec),ncol=p)
    B[upper.tri(B)] <- 0
    d<-d+rho_d*grad_d
    
    # calculate lowerbound
    
    
    # tlowerbound[i]=system.time(Lbound[i]<-cal_Lbound(theta_u=theta,nu=nu,mu=mu,d=d,B=B,m=m,
    # x=x,yo=yo,wpluswt=wpluswt,wtw=wtw))[[3]]
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    
    
    all_paras[i,]<-as.vector(mu)
    # Lbound[i]<-0
    ##############################
    # print(i)
  }
  
  return(list(mu=mu,B=B,d=d,yunote=yunote,
              all_paras=all_paras,accepts=accepts,
              dis_yug_yo=mysamples[[3]]))
} 


MNAR_VB_thetaAug_fullB_updating<-function(x,xm,m,w,yo,startyu,p,   #HVB-AllB algorithm
                                          start_theta,N,N1,bsize){ 
  # p=10
  w<-as(w,"CsparseMatrix")
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=ncol(xm)+1     #number of logistic covariates. i.e. betas
  totpara<-rc+lc+2
  adapt_epsilon=10^(-6)
  v=0.95
  
  n<-length(m)
  no<-length(yo)
  nu<-n-no
  
  
  # bsize=10
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  
  grad_2<-function(theta_u,wpluswt,wtw,w,x,I,yo,xo,xu,xm,m,mu_,nu,no,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p means the dimension of X.
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)] #transform of sigma2
    lamda<-theta_u[(nu+p+2)] #transform of sigma rho
    psi<-theta_u[(nu+p+3):length(theta_u)]
    
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    omega_psi=10000
    
    
    y=c(yo,yu)
    # A=I-rho*w
    # M=t(A)%*%A
    M=I-rho*wpluswt+rho*rho*wtw
    #betas
    div_h_beta<-exp(-gama)*t(y-x%*%beta)%*%M%*%x-t(beta)/omega_beta
    
    #gamma
    div_h_gamma<--0.5*n+(0.5*exp(-gama))*t(y-x%*%beta)%*%M%*%(y-x%*%beta)-gama/omega_gamma
    
    #Lamda
    drho_by_dlamda<-2*exp(lamda)/(1+exp(lamda))^2
    
    dM_by_dlamda<-(-wpluswt+2*rho*wtw)*drho_by_dlamda
    # dlogdetM_by_dlamda<-sum(diag(solve(M)%*%dM_by_dlamda))
    
    # tic()
    # M.in.dMbdlamda<-solve(M,dM_by_dlamda) #1
    # print(paste("direct", sum(diag(M.in.dMbdlamda))))
    # toc()
    
    M.in.dMbdlamda<-solve(Cholesky(M),dM_by_dlamda) #2, faster
    # print(paste("Choleskey", sum(diag(M.in.dMbdlamda))))
    
    dlogdetM_by_dlamda<-sum(diag(M.in.dMbdlamda))
    
    div_h_lamda<-0.5*dlogdetM_by_dlamda-0.5*exp(-gama)*t(y-x%*%beta)%*%dM_by_dlamda%*%(y-x%*%beta)-lamda/omega_lamda
    #Psi
    z<-cbind(xm,y)
    div_h_psi<-(t(m-(exp(z%*%psi)/((1+exp(z%*%psi)))))%*%z)-t(psi)/omega_psi
    
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_gamma),as.vector(div_h_lamda),as.vector(div_h_psi))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  
  make.w_x_<-function(w,x,xm,m,k,reminder){
    
    if(reminder==0){
      
      w.list<-list()
      x.list<-list()
      xm.list<-list()
      m.list<-list()
      
      all.index<-1:n
      o<-1:no
      u<-setdiff(all.index,o)
      current<-0
      for(i in 1:k){ #i=1
        ui<-(no+1+current):(no+current+bsize)
        ui_<-setdiff(u,ui)
        
        w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
        # w[1:15,1:15]
        xui<-matrix(x[ui,],nrow=bsize,)
        xui_<-as.matrix(x[ui_,],nrow=nu-bsize)
        xo<-as.matrix(x[o,])
        x.list[[i]]<-rbind(xui,rbind(xui_,xo))
        
        xmui<-matrix(xm[ui,],nrow=bsize,)
        xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
        xmo<-as.matrix(xm[o,])
        xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
        
        m.list[[i]]<-rev(m)
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        w.list<-list()
        x.list<-list()
        xm.list<-list()
        m.list<-list()
        
        all.index<-1:n
        o<-1:no
        u<-setdiff(all.index,o)
        current<-0
        # i=2 current=5
        for(i in 1:(k+1)){
          if(i<=k){
            ui<-(no+1+current):(no+current+bsize)
            ui_<-setdiff(u,ui)
            
            w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
            # w[1:15,1:15]
            xui<-matrix(x[ui,],nrow=bsize,)
            xui_<-matrix(x[ui_,],nrow=nu-bsize)
            xo<-x[o,]
            x.list[[i]]<-rbind(xui,rbind(xui_,xo))
            
            xmui<-matrix(xm[ui,],nrow=bsize,)
            xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
            xmo<-as.matrix(xm[o,])
            xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
            
            
            m.list[[i]]<-rev(m)
            
            current<-current+bsize
          }else{
            ui<-(no+1+current):(no+current+reminder)
            ui_<-setdiff(u,ui)
            
            w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
            # w[1:15,1:15]
            xui<-matrix(x[ui,],nrow=reminder,)
            xui_<-as.matrix(x[ui_,],nrow=nu-reminder)
            xo<-as.matrix(x[o,])
            x.list[[i]]<-rbind(xui,rbind(xui_,xo))
            
            xmui<-matrix(xm[ui,],nrow=reminder,)
            xmui_<-as.matrix(xm[ui_,],nrow=nu-reminder)
            xmo<-as.matrix(xm[o,])
            xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
            
            
            m.list[[i]]<-rev(m)
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(w.list=w.list,x.list=x.list,xm.list=xm.list,m.list=m.list))
  }
  
  sim_yugyotheta_<-function(yo,listw_x_m,    # After adding un-balanced # of blocks 
                            yu.start,beta,
                            sigma2y,rho,psi,
                            N1,k){
    nu=length(yu.start)
    
    
    cal.con.dis<-function(x,yo,yui_,w,I,beta,rho,sigma2y,bsize,no,nu,reminder,j,k){ #j is the iteration number
      
      M<-t(I-rho*w)%*%(I-rho*w)
      # Sigma<-sigma2y*solve(t(I-rho*w)%*%(I-rho*w))
      if(reminder==0){
        # Sigma_ui.ui_o<-Sigma[1:bsize,(bsize+1):n] #between cov(yui and yui_,yo)
        # Sigma_ui_o.ui_o<-Sigma[(bsize+1):n,(bsize+1):n] #cov(yui_,yo and yui_,yo)
        # Sigma_uiui<-Sigma[1:bsize,1:bsize]
        # 
        # 
        # mu_uigui_o<-x[1:bsize,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
        #                                                      c(yui_,yo)-(x[(bsize+1):n,])%*%beta)
        # sigma_uigui_o<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
        
        
        # sigma_uj.v<-Sigma[1:bsize,(bsize+1):n]
        # sigma_v.v<-Sigma[(bsize+1):n,(bsize+1):n] 
        # sigma_uj.uj<-Sigma[1:bsize,1:bsize] 
        # mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        # v<-c(yui_,yo)
        # 
        # mu_ujguv<-x[1:bsize,]%*%beta+sigma_uj.v%*%solve(sigma_v.v,(v-mu_v))
        # sigma_ujgv<-sigma_uj.uj-sigma_uj.v%*%solve(sigma_v.v,t(sigma_uj.v))
        #######################################
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(yui_,yo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        
        
        # j=3
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the falf block
          # as.matrix(Sigma)[1:reminder,(reminder+1):n,drop = FALSE]
          # Sigma_ui.ui_o<-Sigma[1:reminder,(reminder+1):n,drop = FALSE] #between cov(yui and yui_,yo)
          # Sigma_ui_o.ui_o<-Sigma[(reminder+1):n,(reminder+1):n] #cov(yui_,yo and yui_,yo)
          # Sigma_uiui<-Sigma[1:reminder,1:reminder] 
          # 
          # mu_ujguv<-x[1:reminder,]%*%beta+Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,
          #                                                       c(yui_,yo)-(x[-(1:reminder),])%*%beta)
          # sigma_ujgv<-Sigma_uiui-Sigma_ui.ui_o%*%solve(Sigma_ui_o.ui_o,t(Sigma_ui.ui_o))
          
          ################################
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          
          
        }
        
      }
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yu,yo)
      z=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-z))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    
    yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
    yu.chain[1,]<-yu.start
    current_state<-yu.start
    
    # tic()
    accepts<-0
    for(i in 2:N1){#i=3 N=5
      
      
      # j=1
      current<-0
      for (j in 1:(k+ifelse(reminder>0,1,0))) {
        # accepts<-0
        # yubj<-current_state[(current+1):(current+bsize)]
        if(j<=k){ # for full blocks j=2
          yubj_<-current_state[-((current+1):(current+bsize))]
          
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
          proposal_state<-current_state
          
          proposal_state[(current+1):(current+bsize)]<-yubj.star
          
          # current<-5
          ration=min(0,log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j]],psi=psi)-
                       log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j]],psi=psi))
          
          
          
          if ((!is.na(ration))){
            if((log(runif(1)) <= ration))
            {
              current_state <- proposal_state
              accepts<-accepts+1
            }
            
          }
          
          current<-current+bsize
        }else{ # for the half block j=2 current=5
          yubj_<-current_state[-((current+1):(current+reminder))]
          # current=10
          # print("in half block")
          co.dis<-cal.con.dis(x=listw_x_m$x.list[[j]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j]],I=I,beta=beta,
                              rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j,k=k)
          # print("out half block")
          con.mean<-co.dis$mu_ujguv
          con.var<-co.dis$sigma_ujgv
          
          yubj.star=rmvn(1,as.vector(con.mean),as.matrix((con.var)))
          proposal_state<-current_state
          
          proposal_state[(current+1):(current+reminder)]<-yubj.star
          
          # current<-5
          ration=min(0,log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j]],psi=psi)-
                       log.target(m=listw_x_m$m.list[[j]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j]],psi=psi))
          
          
          
          if ((!is.na(ration))){
            if((log(runif(1)) <= ration))
            {
              current_state <- proposal_state
              accepts<-accepts+1
            }
            
          }
          
        }
      }
      
      
      
      
      yu.chain[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  
  #to calculate lower bound
  #initial values for lamdha
  #p=2
  #mu<-rep(0.1,totpara)
  
  mu<-c(start_theta)
  startyu<-yu
  d<-rep(0.1,totpara)
  B_vec<-matrix(rep(0.1,(totpara)*p),ncol=1)
  B<-matrix(as.vector(B_vec),ncol=p)
  B[upper.tri(B)] <- 0
  # B<-as(B, "dgCMatrix")
  
  #initial values for adaptive learning
  E_g2_mu<-rep(0,totpara)
  E_g2_d<-rep(0,totpara)
  E_g2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  E_delta2_mu<-rep(0,totpara)
  E_delta2_d<-rep(0,totpara)
  E_delta2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  rho_mu<-rep(0,totpara)
  rho_B<-rep(0,totpara*p)
  rho_d<-rep(0,totpara)
  
  # creating a vector for storing lowerbounds
  Lbound<-c()
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) #To store lamda in each iteration
  tlowerbound<-c()
  accepts<-c()
  
  
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder =reminder)
  
  #i=1
  for(i in 1:N){
    
    # thetas
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    z<-matrix(zeta[1:p],ncol = 1)
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    theta<-mu+B%*%z+d*as.vector(epsilon)
    D<-diag(d)
    
    # generate yus from p(yu | yo, theta)
    beta<-theta[(1:rc)]
    gamma.lamda<-theta[(rc+1):(rc+2)]
    rho<-(exp(gamma.lamda[2])-1)/(exp(gamma.lamda[2])+1)
    sigma2y<-exp(gamma.lamda[1])
    psi<-theta[-(1:(rc+2))]
    # N1=100
    
    mysamples<-sim_yugyotheta_(yo=yo,listw_x_m=listw_x_m,
                               yu.start=startyu,beta=beta,
                               sigma2y=sigma2y,rho=rho,
                               psi=psi,N1=N1,k=k)
    # plot(ts((mysamples$samples)[,1]))
    accepts[i]<-mysamples$accepts
    # dim(mysamples$samples)
    # yunote<-apply((mysamples$samples)[-(1:N1/2),], 2, mean)
    yunote<-(mysamples$yugyom.chain)[N1,]
    startyu<-yunote # set starting yu of the MCMC samples of the next iteration
    
    theta_u<-c(as.vector(yunote),as.vector(theta)) # combine thetas and yu
    
    # tic()
    gradg_theta<-grad_2(theta_u=theta_u,w=w,x,I=I,yo=yo,xo=xo,xu=xu,xm=xm,
                        wpluswt=wpluswt,wtw=wtw,nu=nu,no=no,n=n,mu_=mu_,m=m) # gradient of log h(theta)
    
    # toc()
    
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    
    #Construct UB estimates
    grad_mu<-gradg_theta-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradg_theta))
    grad_B<-t(t(z)%x%I_theta)%*%(gradg_theta-dlogq_by_dtheta)
    gradd<-epsilon*(gradg_theta-dlogq_by_dtheta)
    grad_d<-matrix(gradd,ncol = 1)
    
    
    #Set new learning rates
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B_vec<-B_vec+rho_B*grad_B
    B<-matrix(as.vector(B_vec),ncol=p)
    B[upper.tri(B)] <- 0
    d<-d+rho_d*grad_d
    
    # calculate lowerbound
    
    
    # tlowerbound[i]=system.time(Lbound[i]<-cal_Lbound(theta_u=theta,nu=nu,mu=mu,d=d,B=B,m=m,
    # x=x,yo=yo,wpluswt=wpluswt,wtw=wtw))[[3]]
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    
    
    all_paras[i,]<-as.vector(mu)
    # Lbound[i]<-0
    ##############################
    # print(i)
  }
  
  return(list(mu=mu,B=B,d=d,yunote=yunote,
              all_paras=all_paras,accepts=accepts,
              dis_yug_yo=mysamples[[3]]))
} 




MNAR_VB_thetaAug_R3B_updating<-function(x,xm,m,w,yo,startyu,p,    #HVB-3B algorithm
                                        start_theta,N,N1,bsize){ 
  # p=10
  w<-as(w,"CsparseMatrix")
  x<-cbind(1,x)
  xm<-cbind(1,xm)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=ncol(xm)+1     #number of logistic covariates. i.e. betas
  totpara<-rc+lc+2
  adapt_epsilon=10^(-6)
  v=0.95
  
  n<-length(m)
  no<-length(yo)
  nu<-n-no
  
  
  
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  I<-Diagonal(n)
  I_p<-Diagonal(p)
  
  
  xo=x[1:no,]
  xu=x[-(1:no),]
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  
  grad_2<-function(theta_u,wpluswt,wtw,w,x,I,yo,xo,xu,xm,m,mu_,nu,no,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p means the dimension of X.
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)] #transform of sigma2
    lamda<-theta_u[(nu+p+2)] #transform of sigma rho
    psi<-theta_u[(nu+p+3):length(theta_u)]
    
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    omega_psi=10000
    
    
    y=c(yo,yu)
    # A=I-rho*w
    # M=t(A)%*%A
    M=I-rho*wpluswt+rho*rho*wtw
    #betas
    div_h_beta<-exp(-gama)*t(y-x%*%beta)%*%M%*%x-t(beta)/omega_beta
    
    #gamma
    div_h_gamma<--0.5*n+(0.5*exp(-gama))*t(y-x%*%beta)%*%M%*%(y-x%*%beta)-gama/omega_gamma
    
    #Lamda
    drho_by_dlamda<-2*exp(lamda)/(1+exp(lamda))^2
    
    dM_by_dlamda<-(-wpluswt+2*rho*wtw)*drho_by_dlamda
    # dlogdetM_by_dlamda<-sum(diag(solve(M)%*%dM_by_dlamda))
    
    # tic()
    # M.in.dMbdlamda<-solve(M,dM_by_dlamda) #1
    # print(paste("direct", sum(diag(M.in.dMbdlamda))))
    # toc()
    
    M.in.dMbdlamda<-solve(Cholesky(M),dM_by_dlamda) #2, faster
    # print(paste("Choleskey", sum(diag(M.in.dMbdlamda))))
    
    dlogdetM_by_dlamda<-sum(diag(M.in.dMbdlamda))
    
    div_h_lamda<-0.5*dlogdetM_by_dlamda-0.5*exp(-gama)*t(y-x%*%beta)%*%dM_by_dlamda%*%(y-x%*%beta)-lamda/omega_lamda
    #Psi
    z<-cbind(xm,y)
    div_h_psi<-(t(m-(exp(z%*%psi)/((1+exp(z%*%psi)))))%*%z)-t(psi)/omega_psi
    
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_gamma),as.vector(div_h_lamda),as.vector(div_h_psi))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  
  make.w_x_<-function(w,x,xm,m,k,reminder){
    
    if(reminder==0){
      
      w.list<-list()
      x.list<-list()
      xm.list<-list()
      m.list<-list()
      
      all.index<-1:n
      o<-1:no
      u<-setdiff(all.index,o)
      current<-0
      for(i in 1:k){ #i=1
        ui<-(no+1+current):(no+current+bsize)
        ui_<-setdiff(u,ui)
        
        w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
        # w[1:15,1:15]
        xui<-matrix(x[ui,],nrow=bsize,)
        xui_<-as.matrix(x[ui_,],nrow=nu-bsize)
        xo<-as.matrix(x[o,])
        x.list[[i]]<-rbind(xui,rbind(xui_,xo))
        
        xmui<-matrix(xm[ui,],nrow=bsize,)
        xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
        xmo<-as.matrix(xm[o,])
        xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
        
        m.list[[i]]<-rev(m)
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        w.list<-list()
        x.list<-list()
        xm.list<-list()
        m.list<-list()
        
        all.index<-1:n
        o<-1:no
        u<-setdiff(all.index,o)
        current<-0
        # i=2 current=5
        for(i in 1:(k+1)){
          if(i<=k){
            ui<-(no+1+current):(no+current+bsize)
            ui_<-setdiff(u,ui)
            
            w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
            # w[1:15,1:15]
            xui<-matrix(x[ui,],nrow=bsize,)
            xui_<-matrix(x[ui_,],nrow=nu-bsize)
            xo<-x[o,]
            x.list[[i]]<-rbind(xui,rbind(xui_,xo))
            
            xmui<-matrix(xm[ui,],nrow=bsize,)
            xmui_<-as.matrix(xm[ui_,],nrow=nu-bsize)
            xmo<-as.matrix(xm[o,])
            xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
            
            
            m.list[[i]]<-rev(m)
            
            current<-current+bsize
          }else{
            ui<-(no+1+current):(no+current+reminder)
            ui_<-setdiff(u,ui)
            
            w.list[[i]]<-w[c(ui,ui_,o),c(ui,ui_,o)]
            # w[1:15,1:15]
            xui<-matrix(x[ui,],nrow=reminder,)
            xui_<-as.matrix(x[ui_,],nrow=nu-reminder)
            xo<-as.matrix(x[o,])
            x.list[[i]]<-rbind(xui,rbind(xui_,xo))
            
            xmui<-matrix(xm[ui,],nrow=reminder,)
            xmui_<-as.matrix(xm[ui_,],nrow=nu-reminder)
            xmo<-as.matrix(xm[o,])
            xm.list[[i]]<-rbind(xmui,rbind(xmui_,xmo))
            
            
            m.list[[i]]<-rev(m)
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(w.list=w.list,x.list=x.list,xm.list=xm.list,m.list=m.list))
  }
  
  
  sim_yugyotheta_<-function(yo,make.w_x_,    # After adding un-balanced # of blocks 
                            yu.start,beta,
                            sigma2y,rho,psi,
                            N1,k,bsize,reminder){
    nu=length(yu.start)
    
    
    cal.con.dis<-function(x,yo,yui_,w,I,beta,rho,sigma2y,bsize,no,nu,reminder,j,k){ #j is the iteration number
      
      M<-t(I-rho*w)%*%(I-rho*w)
      
      if(reminder==0){
        M_v.v=M[(bsize+1):n,(bsize+1):n]
        M_uj.uj=M[1:bsize,1:bsize]
        M_uj.v=M[1:bsize,(bsize+1):n]
        M_v.uj=t(M_uj.v)
        v<-c(yui_,yo)
        mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
        
        In<-Diagonal(ncol(M_uj.uj))
        cholMuj.uj<-Cholesky(M_uj.uj)
        mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
        sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
        
      }else{
        if(j<=k){
          M_v.v=M[(bsize+1):n,(bsize+1):n]
          M_uj.uj=M[1:bsize,1:bsize]
          M_uj.v=M[1:bsize,(bsize+1):n]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c(x[(bsize+1):nu,]%*%beta,x[-(1:nu),]%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:bsize,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          # dim(t(Sigma_ui.ui_o))
          
        }else{ # for the falf block
          
          M_v.v=M[(reminder+1):n,(reminder+1):n]
          M_uj.uj=M[1:reminder,1:reminder]
          M_uj.v=M[1:reminder,(reminder+1):n,drop = FALSE]
          M_v.uj=t(M_uj.v)
          v<-c(yui_,yo)
          mu_v<-c((x[-(1:reminder),])%*%beta)
          
          In<-Diagonal(ncol(M_uj.uj))
          cholMuj.uj<-Cholesky(M_uj.uj)
          mu_ujguv=x[1:reminder,]%*%beta-solve(cholMuj.uj,M_uj.v)%*%(v-mu_v)
          sigma_ujgv=forceSymmetric(sigma2y*solve(cholMuj.uj,In))
          
        }
        
      }
      return(list(mu_ujguv=mu_ujguv,sigma_ujgv=sigma_ujgv))
    }
    log.target<-function(m,yo,yu,x,psi){
      
      y=c(yu,yo)
      z=cbind(x,y)%*%psi
      pxy<-1/(1+exp(-z))
      
      loglike1<-sum(m*log(pxy)+(1-m)*log(1-pxy)) #p_m
      return(loglike1)
    }
    
    j.list<-sample(1:length(listw_x_m$x.list),size = 3,replace =F) # selecting a block to be updated
    accepts<-0
    if(length(which(j.list ==  k+1)) >0){ # check wether, in the jlist, the half block in included 
      for(i in 1:N1){
        yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
        yu.chain[1,]<-yu.start
        current_state<-yu.start
        
        yubj_<-current_state[(1:(k*bsize))]
        co.dis<-cal.con.dis(x=listw_x_m$x.list[[k+1]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[k+1]],I=I,beta=beta,
                            rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=k+1,k=k)
        
        con.mean.halfj<-co.dis$mu_ujguv
        con.var.halfj<-co.dis$sigma_ujgv
        
        yubj.star=rmvn(1,as.vector(con.mean.halfj),as.matrix((con.var.halfj)))
        proposal_state<-current_state
        proposal_state[(k*bsize+1):(k*bsize+reminder)]<-yubj.star
        
        ration=min(0,log.target(m=listw_x_m$m.list[[k+1]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[k+1]],psi=psi)-
                     log.target(m=listw_x_m$m.list[[k+1]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[k+1]],psi=psi))
        
        
        if ((!is.na(ration))){
          if((log(runif(1)) <= ration))
          {
            current_state <- proposal_state
            accepts<-accepts+1
          }
          
        }
        
        j.list<- j.list[j.list !=  k+1] # selecting block numbers other than the half block
        
        
        # random full block1
        rj1.start<-bsize*(j.list[1]-1)+1 #starting index of this block
        
        yubj_<-current_state[-(rj1.start:(rj1.start-1+bsize))]
        
        co.dis<-cal.con.dis(x=listw_x_m$x.list[[j.list[1]]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j.list[1]]],I=I,beta=beta,
                            rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j.list[1],k=k)
        
        con.mean.jr1<-co.dis$mu_ujguv
        con.var.jr1<-co.dis$sigma_ujgv
        
        
        yubj.star=rmvn(1,as.vector(con.mean.jr1),as.matrix((con.var.jr1)))
        proposal_state<-current_state
        
        proposal_state[rj1.start:(rj1.start-1+bsize)]<-yubj.star
        
        
        ration=min(0,log.target(m=listw_x_m$m.list[[j.list[1]]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j.list[1]]],psi=psi)-
                     log.target(m=listw_x_m$m.list[[j.list[1]]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j.list[1]]],psi=psi))
        
        
        
        if ((!is.na(ration))){
          if((log(runif(1)) <= ration))
          {
            current_state <- proposal_state
            accepts<-accepts+1
          }
          
        }
        
        
        # random full block2
        rj2.start<-bsize*(j.list[2]-1)+1 #starting index of this block
        
        yubj_<-current_state[-(rj2.start:(rj2.start-1+bsize))]
        
        co.dis<-cal.con.dis(x=listw_x_m$x.list[[j.list[2]]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j.list[2]]],I=I,beta=beta,
                            rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j.list[2],k=k)
        
        con.mean.jr2<-co.dis$mu_ujguv
        con.var.jr2<-co.dis$sigma_ujgv
        
        
        yubj.star=rmvn(1,as.vector(con.mean.jr2),as.matrix((con.var.jr2)))
        proposal_state<-current_state
        
        proposal_state[rj2.start:(rj2.start-1+bsize)]<-yubj.star
        
        
        ration=min(0,log.target(m=listw_x_m$m.list[[j.list[2]]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j.list[2]]],psi=psi)-
                     log.target(m=listw_x_m$m.list[[j.list[2]]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j.list[2]]],psi=psi))
        
        
        
        if ((!is.na(ration))){
          if((log(runif(1)) <= ration))
          {
            current_state <- proposal_state
            accepts<-accepts+1
          }
          
        }
        
        yu.chain[i,] <- current_state
        
      }
    }else{
      
      for(i in 1:N1){
        yu.chain<-matrix(rep(0,nu*N1),ncol=nu)
        yu.chain[1,]<-yu.start
        current_state<-yu.start
        
        # random block 1
        rj1.start<-bsize*(j.list[1]-1)+1 #starting index of this block
        
        yubj_<-current_state[-(rj1.start:(rj1.start-1+bsize))]
        
        co.dis<-cal.con.dis(x=listw_x_m$x.list[[j.list[1]]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j.list[1]]],I=I,beta=beta,
                            rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j.list[1],k=k)
        
        con.mean.jr1<-co.dis$mu_ujguv
        con.var.jr1<-co.dis$sigma_ujgv
        
        
        yubj.star=rmvn(1,as.vector(con.mean.jr1),as.matrix((con.var.jr1)))
        proposal_state<-current_state
        
        proposal_state[rj1.start:(rj1.start-1+bsize)]<-yubj.star
        
        
        ration=min(0,log.target(m=listw_x_m$m.list[[j.list[1]]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j.list[1]]],psi=psi)-
                     log.target(m=listw_x_m$m.list[[j.list[1]]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j.list[1]]],psi=psi))
        
        
        
        if ((!is.na(ration))){
          if((log(runif(1)) <= ration))
          {
            current_state <- proposal_state
            accepts<-accepts+1
          }
          
        }
        
        # random block 2
        
        rj1.start<-bsize*(j.list[2]-1)+1 #starting index of this block
        
        yubj_<-current_state[-(rj1.start:(rj1.start-1+bsize))]
        
        co.dis<-cal.con.dis(x=listw_x_m$x.list[[j.list[2]]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j.list[2]]],I=I,beta=beta,
                            rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j.list[2],k=k)
        
        con.mean.jr1<-co.dis$mu_ujguv
        con.var.jr1<-co.dis$sigma_ujgv
        
        
        yubj.star=rmvn(1,as.vector(con.mean.jr1),as.matrix((con.var.jr1)))
        proposal_state<-current_state
        
        proposal_state[rj1.start:(rj1.start-1+bsize)]<-yubj.star
        
        
        ration=min(0,log.target(m=listw_x_m$m.list[[j.list[2]]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j.list[2]]],psi=psi)-
                     log.target(m=listw_x_m$m.list[[j.list[2]]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j.list[2]]],psi=psi))
        
        
        
        if ((!is.na(ration))){
          if((log(runif(1)) <= ration))
          {
            current_state <- proposal_state
            accepts<-accepts+1
          }
          
        }
        
        
        # random full block3
        
        rj2.start<-bsize*(j.list[3]-1)+1 #starting index of this block
        
        yubj_<-current_state[-(rj2.start:(rj2.start-1+bsize))]
        
        co.dis<-cal.con.dis(x=listw_x_m$x.list[[j.list[3]]],yo=yo,yui_=yubj_,w=listw_x_m$w.list[[j.list[3]]],I=I,beta=beta,
                            rho=rho,sigma2y=sigma2y,bsize=bsize,no=no,nu=nu,reminder=reminder,j=j.list[3],k=k)
        
        con.mean.jr2<-co.dis$mu_ujguv
        con.var.jr2<-co.dis$sigma_ujgv
        
        
        yubj.star=rmvn(1,as.vector(con.mean.jr2),as.matrix((con.var.jr2)))
        proposal_state<-current_state
        
        proposal_state[rj2.start:(rj2.start-1+bsize)]<-yubj.star
        
        
        ration=min(0,log.target(m=listw_x_m$m.list[[j.list[3]]],yo=yo,yu=proposal_state,x=listw_x_m$xm.list[[j.list[3]]],psi=psi)-
                     log.target(m=listw_x_m$m.list[[j.list[3]]],yo=yo,yu=current_state,x=listw_x_m$xm.list[[j.list[3]]],psi=psi))
        
        
        
        if ((!is.na(ration))){
          if((log(runif(1)) <= ration))
          {
            current_state <- proposal_state
            accepts<-accepts+1
          }
          
        }
        
        
        
        
        yu.chain[i,] <- current_state
        
      }
      
      
    }
    
    
    updated.block.no<-0
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis,updated.block.no=updated.block.no))
    
  }
  
  
  mu<-c(start_theta)
  startyu<-yu
  d<-rep(0.1,totpara)
  B_vec<-matrix(rep(0.1,(totpara)*p),ncol=1)
  B<-matrix(as.vector(B_vec),ncol=p)
  B[upper.tri(B)] <- 0
  # B<-as(B, "dgCMatrix")
  
  #initial values for adaptive learning
  E_g2_mu<-rep(0,totpara)
  E_g2_d<-rep(0,totpara)
  E_g2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  E_delta2_mu<-rep(0,totpara)
  E_delta2_d<-rep(0,totpara)
  E_delta2_B<-rep(0,totpara*p) # this needs to be changed when p>1
  
  rho_mu<-rep(0,totpara)
  rho_B<-rep(0,totpara*p)
  rho_d<-rep(0,totpara)
  
  # creating a vector for storing lowerbounds
  Lbound<-c()
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) #To store lamda in each iteration
  tlowerbound<-c()
  accepts<-c()
  updated.block.no<-c()
  
  listw_x_m<-make.w_x_(w=w,x=x,xm=xm,m=m,k=k,reminder =reminder)
  
  #i=1
  for(i in 1:N){
    
    # thetas
    zeta<-(rnorm(n=(p+totpara)))
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    z<-matrix(zeta[1:p],ncol = 1)
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    theta<-mu+B%*%z+d*as.vector(epsilon)
    D<-diag(d)
    
    # generate yus from p(yu | yo, theta)
    beta<-theta[(1:rc)]
    gamma.lamda<-theta[(rc+1):(rc+2)]
    rho<-(exp(gamma.lamda[2])-1)/(exp(gamma.lamda[2])+1)
    sigma2y<-exp(gamma.lamda[1])
    psi<-theta[-(1:(rc+2))]
    # N1=100
    
    mysamples<-sim_yugyotheta_(yo=yo,make.w_x_=make.w_x_,
                               yu.start=startyu,beta=beta,
                               sigma2y=sigma2y,rho=rho,
                               psi=psi,N1=N1,k=k,bsize=bsize,reminder=reminder)
    # plot(ts((mysamples$samples)[,1]))
    accepts[i]<-mysamples$accepts
    # dim(mysamples$samples)
    # yunote<-apply((mysamples$samples)[-(1:N1/2),], 2, mean)
    yunote<-(mysamples$yugyom.chain)[N1,]
    startyu<-yunote # set starting yu of the MCMC samples of the next iteration
    updated.block.no[i]<-mysamples$updated.block.no
    theta_u<-c(as.vector(yunote),as.vector(theta)) # combine thetas and yu
    
    # tic()
    gradg_theta<-grad_2(theta_u=theta_u,w=w,x,I=I,yo=yo,xo=xo,xu=xu,xm=xm,
                        wpluswt=wpluswt,wtw=wtw,nu=nu,no=no,n=n,mu_=mu_,m=m) # gradient of log h(theta)
    
    # toc()
    
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    
    #Construct UB estimates
    grad_mu<-gradg_theta-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradg_theta))
    grad_B<-t(t(z)%x%I_theta)%*%(gradg_theta-dlogq_by_dtheta)
    gradd<-epsilon*(gradg_theta-dlogq_by_dtheta)
    grad_d<-matrix(gradd,ncol = 1)
    
    
    #Set new learning rates
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B_vec<-B_vec+rho_B*grad_B
    B<-matrix(as.vector(B_vec),ncol=p)
    B[upper.tri(B)] <- 0
    d<-d+rho_d*grad_d
    
    # calculate lowerbound
    
    
    # tlowerbound[i]=system.time(Lbound[i]<-cal_Lbound(theta_u=theta,nu=nu,mu=mu,d=d,B=B,m=m,
    # x=x,yo=yo,wpluswt=wpluswt,wtw=wtw))[[3]]
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    
    
    all_paras[i,]<-as.vector(mu)
    # Lbound[i]<-0
    ##############################
    # print(i)
  }
  
  return(list(mu=mu,B=B,d=d,yunote=yunote,
              all_paras=all_paras,accepts=accepts,
              dis_yug_yo=mysamples[[3]],updated.block.no=updated.block.no))
} 




