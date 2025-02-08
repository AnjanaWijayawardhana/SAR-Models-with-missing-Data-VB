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
library("vctrs")
library(dplyr)
library(tidyr)
library(reshape2)



######################### Supporting Functions

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




splitter<-function(x,y,w,p_sample){
  
  n=nrow(w)
  ns=round(n*p_sample)
  nu=n-ns
  x=as.matrix(x,nrow=n)
  
  s<-sample(1:n,ns) #sampled
  s<-sort(s)
  u<-setdiff(1:n,s)  # non-sample units
  
  ns<-length(s)
  nu<-n-ns
  
  ys=y[s] #Sold unsold y
  yu=y[u]
  
  xs=as.matrix(x[s,],nrow=ns) #sold unsold x
  xu=as.matrix(x[u,],nrow=nu)
  
  
  w<-w[c(s,u),c(s,u)] # Create w partitioned matrix
  y<-cbind(c(ys,yu)) # create full (observed+unobserved) x
  x<-xs
  x<-rbind(x,xu) # create full (observed+unobserved) y
  
  l=list("W"=w,"Y"=y,"X"=x,"ns"=ns)
  return(l)
  
}   #splitting data set randomly
#splitting data set randomly recommend for SEM


# Functions to find starting values

estimate_regression<-function(x,yo,no){
  # x<-cbind(1,x)
  fit.regression<-lm(yo~x[1:no,])
  return(fit.regression)
}

sim_yugyo<-function(yo,x,w,theta){
  x<-cbind(1,x)
  n<-ncol(w)
  no<-length(yo)
  missingindexes<-(no+1):n
  xu<-x[missingindexes,]
  xo<-x[setdiff(1:n,missingindexes),]
  # x<-ncol(x)
  beta<-paras.st[1:ncol(x)]
  sigma2y<-exp(paras.st[ncol(x)+1])
  rho<-(exp(paras.st[ncol(x)+2]) -1)/(exp(paras.st[ncol(x)+2]) +1)
  
  I<-Diagonal(n)
  I_nu<-Diagonal(n-no)
  A<-I-rho*w
  M<-t(A)%*%A
  M_oo=M[1:no,1:no]
  M_uu=M[(no+1):n,(no+1):n]
  M_ou=M[1:no,(no+1):n]
  M_uo=t(M_ou)
  
  cholMuu<-Cholesky(M_uu)
  mu_ugo<-xu%*%beta-solve(cholMuu,M_uo)%*%(yo-xo%*%beta) 
  sigma_ugo<-sigma2y*solve(cholMuu,I_nu) 
  yu.star<-rmvn(1,as.vector(mu_ugo),as.matrix(sigma_ugo)) 
  return(yu.star)
}


# Functions for sampling from posteriors


simulate.VB.MAR_vb1<-function(mu,B,x,d,p,nu){
  
  
  
  m<-nu+ncol(x)+1+2
  
  epsilon<-rnorm(m)
  z<-rnorm(p)
  
  theta<-mu+B%*%z+as.vector(d)*epsilon
  
  #rho
  lamda<-theta[m]
  rho<-(exp(lamda) -
          1)/(exp(lamda) +1)
  
  #sigma2y
  sigma2y<-exp(theta[m-1])
  
  theta[m]<-rho
  theta[m-1]<-sigma2y
  return(as.vector(theta))
  
} #  JVB


simulate.VB.SEM_vb2<-function(mu,B,x,d,p){
  
  
  nu=0
  m<-ncol(x)+1+2 # missing vals, betas(col.x+1), rho+sigma2y+  psis
  
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
  
  # yu<-theta_u[1:nu]
  # beta<-theta_u[(nu+1):(nu+p)]
  # gama<-theta_u[(nu+p+1)]
  # lamda<-theta_u[(nu+p+2)]
  # psi<-theta_u[(nu+p+2+1):length(theta_u)]
  # rho<-(exp(lamda)-1)/(exp(lamda)+1)
  
  
  
} # HVB 





generate_post_yu_MCMC_B<-function(theta,w.split,x.split,yo,bsize,yu0,N1){
  # 1 calculate det(M) usigng chol(), and 2 uses choleskey(). 2 is very fast
  make.w_x_<-function(w,x,k,reminder){
    
    if(reminder==0){
      
      w.list<-list()
      x.list<-list()
      
      
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
        
        
        
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        w.list<-list()
        x.list<-list()
        
        
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
            
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(w.list=w.list,x.list=x.list))
  }
  
  sim_yugyotheta_MCMC_G<-function(yo,listw_x_m,    # After adding un-balanced # of blocks 
                                  yu.start,theta,
                                  N1,k,bsize,reminder){
    
    beta=theta[1:(ncol(listw_x_m$x.list[[1]]))]
    sigma2y=theta[ncol(listw_x_m$x.list[[1]])+1]
    rho=theta[ncol(listw_x_m$x.list[[1]])+2]
    
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
          
          current_state <- proposal_state
          
          
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
          current_state <- proposal_state
          accepts<-accepts+1
          
          
        }
      }
      
      
      
      
      yu.chain[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  # bsize=500
  nu<-length(yu0)
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  
  x.split<-cbind(1,x.split)
  
  
  n<-nrow(w.split)
  # nu<-sum(is.na(y_ou))
  no<-n-nu
  I<-Diagonal(n)
  missingindexes<-(no+1):n
  xu<-x.split[missingindexes,]
  xo<-x.split[setdiff(1:n,missingindexes),]
  
  
  
  
  listw_x_m<-make.w_x_(w=w.split,x=x.split,k=k,reminder=reminder)
  
  
  
  yu.star.result<-sim_yugyotheta_MCMC_G(yo=yo,listw_x_m=listw_x_m,
                                        yu.start=yu0,theta
                                        ,N1=N1,k=k,bsize,reminder)
  return(yu.star.result$yugyom.chain[N1,])
  
}


######################### VB and MML (Suesse, 2018) algorithms


SEM_MAR_VB.new<-function(x,N,w,y,no,
                         p,start_theta){ 
  # p=1
  x<-cbind(1,x)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=0     #number of logistic covariates. i.e. betas
  totpara<-rc+lc+2
  adapt_epsilon=10^(-6)
  v=0.95
  I_p<-Diagonal(p)
  
  n<-length(y)
  # nu<-sum(is.na(y_ou))
  nu<-n-no
  I<-Diagonal(n)
  
  missingindexes<-(no+1):n
  xu<-x[missingindexes,]
  xo<-x[setdiff(1:n,missingindexes),]
  yo<-y[setdiff(1:n,missingindexes)]
  
  
  # u<-missingindexes
  # o<-setdiff(1:n,u)
  # 
  # w<-w[c(o,u),c(o,u)]
  # x<-rbind(xo,xu)
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  # nu<-3
  # theta_u<-c(rep(0,nu),c(1,2),0.5,0.1,c(3,3,3))
  grad<-function(theta_u,wpluswt,wtw,w,x,xo,xu,I,yo,nu,no,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p means the dimension of X.
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)] #transform of sigma2
    lamda<-theta_u[(nu+p+2)] #transform of sigma rho
    
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    
    
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
    
    #############################- M Direct inverse
    # tic()
    # dlogdetM_by_dlamda<-sum(diag(solve(M)%*%dM_by_dlamda))
    # print(dlogdetM_by_dlamda)
    # toc()
    
    # print(class(M))
    # print(class(dM_by_dlamda))
    #############################- M Choleskey inverse
    # tic()
    dlogdetM_by_dlamda<-sum(diag(solve(Cholesky(M),dM_by_dlamda)))
    # print(dlogdetM_by_dlamda)
    # toc()
    ########################
    
    div_h_lamda<-0.5*dlogdetM_by_dlamda-0.5*exp(-gama)*t(y-x%*%beta)%*%dM_by_dlamda%*%(y-x%*%beta)-lamda/omega_lamda
    #Psi
    M_uu<-M[(no+1):n,(no+1):n]
    M_uo<-M[(no+1):n,1:no]
    # M_ou<-M[1:no,(no+1):n]
    ru<-yu-xu%*%beta
    ro<-yo-xo%*%beta
    div_h_yu_2<--exp(-gama)*(0.5*t(M_uo%*%ro)+t(ru)%*%M_uu)
    
    #3
    
    # M_oo<-M[1:no,1:no]
    # M_uu<-M[(no+1):n,(no+1):n]
    # M_uo<-M[(no+1):n,1:no]
    # M_ou<-M[1:no,(no+1):n]
    # ru<-yu-xu%*%beta
    # ro<-yo-xo%*%beta
    # div_h_yu_2<--0.5*exp(-gama)*(t(M_uo%*%ro)+t(ro)%*%M_ou+2*t(ru)%*%M_uu)
    
    # dim(div_h_yu_2)
    
    div_h_yu<-div_h_yu_2
    
    fullgrad<-c(as.vector(div_h_yu),as.vector(div_h_beta),
                as.vector(div_h_gamma),as.vector(div_h_lamda))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  # nu<-3
  # theta_u<-c(rep(1,3),-1,-2,0.1,0.2,rep(9,3))
  
  cal_Lbound<-function(theta_u,nu,mu,d,B,x,yo,wpluswt,wtw){
    
    p<-ncol(x) #length of betas
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)]
    lamda<-theta_u[(nu+p+2)]
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    y<-c(yo,yu)
    n<-length(y)
    A=I-rho*w
    M=t(A)%*%A
    
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    D<-diag(as.vector(d))
    
    
    
    
    # sigma_<-exp(gama)*solve(M)
    # loglike2<-dmvnorm(x=as.vector(y),mean=x%*%beta,
    #                   sigma=sigma_, log=T,checkSymmetry = F) # direct p_y 
    
    log_det<-2*sum(log(diag(chol(M))))
    loglike2<--0.5*n*log(2*pi)-0.5*n*gama+
      0.5*log_det-0.5*exp(-gama)*t(y-x%*%beta)%*%M%*%(y-x%*%beta) # p_y efficient
    
    loglike3<-dmvnorm(x=as.vector(beta),mean = rep(0,ncol(x)),sigma=omega_beta*diag(ncol(x)),
                      log=T,checkSymmetry = F) #p_beta
    loglike4<-dnorm(x=gama,mean = 0,sd=omega_gamma) #p_gama (variance)
    
    loglike5<-dnorm(x=lamda,mean = 0,sd=omega_lamda) #p_lamda (variance)
    
    
    logq<-dmvnorm(x=as.vector(theta_u),mean = as.vector(mu),sigma =B%*%t(B)+
                    D^2,log = T,checkSymmetry = F) # q
    
    
    return(as.numeric(loglike2)+
             as.numeric(loglike3)+as.numeric(loglike4)+
             as.numeric(loglike5)-as.numeric(logq))
    
    
  } #to calculate lower bound
  
  cal_Lbound_2<-function(theta_u,nu,mu,d,B,x,yo,wpluswt,wtw){
    
    p<-ncol(x) #length of betas
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)]
    lamda<-theta_u[(nu+p+2)]
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    
    p_<-ncol(B)
    y<-c(yo,yu)
    n<-length(y)
    M<-I-rho*wpluswt+rho^2*wtw
    nu_theta<-length(theta_u)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_psi=10000
    omega_lamda=10000
    
    D<-Diagonal(n=length(d),x=as.vector(d)) #making D sparse and diagonal
    
    
    
    # log_det<-2*sum(log(diag(chol(M))))
    
    log_det<-2*determinant(Cholesky(M))[[1]][[1]]
    loglike2<--0.5*n*log(2*pi)-0.5*n*gama+
      0.5*log_det-0.5*exp(-gama)*t(y-x%*%beta)%*%M%*%(y-x%*%beta) # p_y efficient
    
    loglike3<-dmvnorm(x=as.vector(beta),mean = rep(0,ncol(x)),sigma=omega_beta*diag(ncol(x)),
                      log=T,checkSymmetry = F) #p_beta
    loglike4<-dnorm(x=gama,mean = 0,sd=omega_gamma) #p_gama (variance)
    
    loglike5<-dnorm(x=lamda,mean = 0,sd=omega_lamda) #p_lamda (variance)
    
    # logq<-dmvnorm(x=as.vector(theta_u),mean = as.vector(mu),sigma =B%*%t(B)+D^2,log = T,checkSymmetry = F) # q1
    
    ################################################################################ Q using Davids method
    
    r<-theta_u-mu
    B_tilda<-(1/d)*B
    r_tilde<- (1/d)*r
    
    epsilon_hat<-t(B_tilda)%*%r_tilde
    I_BtB<-forceSymmetric(Diagonal(n=p_,x=1)+t(B_tilda)%*%B_tilda) # o/w, i.e if not "forceSymmetric" logq may generate NaN values 
    
    ######## chol
    R<-as(chol(I_BtB),"CsparseMatrix")
    r_tilde2<-solve(t(R),epsilon_hat)
    rtcovr<-t(r_tilde)%*%r_tilde-t(r_tilde2)%*%r_tilde2 #Davids
    logdet<-sum(log(as.vector(d^2)))+sum(log(diag(R^2)))
    
    ####### choleskey
    # R<-(Matrix::Cholesky(I_BtB)) We casnt use Cholesky since I_BtB is not sparse.
    
    
    logq<--0.5*nu_theta*log(2*pi)-0.5*logdet-0.5*rtcovr #q2
    
    
    ##########################################################################
    
    # logdet<-2*sum(log(diag(chol(B%*%t(B)+D^2))))
    # logq<--0.5*nu_theta*log(2*pi)-0.5*logdet-0.5*t(r)%*%solve(B%*%t(B)+D^2,r) #q3
    
    
    
    
    return(as.numeric(loglike2)+
             as.numeric(loglike3)+as.numeric(loglike4)+
             as.numeric(loglike5)-as.numeric(logq))
    
    
  } #to calculate lower bound
  
  # 1 calculate det(M) usigng chol(), and 2 uses choleskey(). 2 is very fast
  
  #initial values for lamdha
  #p=2
  # mu<-rep(0.1,nu)
  mu<-c(start_theta)
  d<-rep(0.1,totpara+nu)
  B<-matrix(rep(0.1,(totpara+nu)*p),ncol=p)
  
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
  
  ###################################### can be removed
  RMS_delta_mu<-rep(0,totpara+nu)
  RMS_delta_B<-rep(0,totpara+nu)
  RMS_delta_d<-rep(0,totpara+nu)
  ####################################
  
  # creating a vector for storing lowerbounds
  Lbound<-c()
  Lbound_2<-c()
  alllamda<-c() #To store lamda in each iteration
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) 
  #i=1
  for(i in 1:N){
    
    zeta<-(rnorm(n=(p+totpara+nu)))
    
    z<-matrix(zeta[1:p],ncol = 1)
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    
    
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    
    
    theta<-mu+B%*%z+d*as.vector(epsilon)
    # theta_t<-matrix(theta_t,ncol=1)
    D<-diag(d)
    # dim(D)
    # Calculate log gradient of theta
    
    gradh_theta<-grad(theta_u=theta,w=w,x,xo=xo,xu=xu,I=I,yo=yo,
                      wpluswt=wpluswt,wtw=wtw,nu=nu,no=no,n=n) # gratient of log h(theta)
    
    # ###-part-2 cal1
    # tic()
    # part2<-solve(forceSymmetric((B%*%t(B)+D^2)))%*%(B%*%z+d*as.vector(epsilon))# the second part inside the expectation
    # toc()                # Direct calculation
    # print(part2[2:4])
    
    ###-part-2 cal2
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+
                                            t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    
    # 
    #Construct UB estimates
    grad_mu<-gradh_theta+part2
    grad_B<-gradh_theta%*%t(z)+part2%*%t(z)
    gradd<-diag(gradh_theta%*%t(epsilon)+part2%*%t(epsilon))
    grad_d<-matrix(gradd,ncol = 1)
    
    #Set new learning rates
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    # RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    # RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    # RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B<-B+rho_B*grad_B
    d<-d+rho_d*grad_d
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    ######################################## can be removed
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    #############################################
    
    # calculate lowerbound
    # tic()
    # Lbound[i]<-cal_Lbound(theta_u=theta,nu=nu,mu=mu,d=d,B=B,
    #                       x=x,yo=yo,wpluswt=wpluswt,wtw=wtw)
    # toc()
    # tic()
    Lbound[i]<-cal_Lbound_2(theta_u=theta,nu=nu,mu=mu,d=d,B=B,
                            x=x,yo=yo,wpluswt=wpluswt,wtw=wtw)
    # toc()
    all_paras[i,]<-as.vector(mu)
    # Lbound[i]<-0
    ##############################
  }
  
  return(list(mu=mu,B=B,d=d,Lbound=Lbound,all_paras=all_paras))
} # JVB


SEM_MAR_VB_aug_ADADELTA<-function(x,N,w,yo,no,
                                  p,start_theta){ 
  # p=1
  x<-cbind(1,x)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=0     #number of logistic covariates. i.e. betas
  totpara<-rc+lc+2
  adapt_epsilon=10^(-6)
  v=0.95
  I_p<-Diagonal(p)
  
  n<-length(y)
  # nu<-sum(is.na(y_ou))
  nu<-n-no
  I<-Diagonal(n)
  I_nu<-Diagonal(nu)
  missingindexes<-(no+1):n
  xu<-x[missingindexes,]
  xo<-x[setdiff(1:n,missingindexes),]
  
  # u<-missingindexes
  # o<-setdiff(1:n,u)
  # 
  # w<-w[c(o,u),c(o,u)]
  # x<-rbind(xo,xu)
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  # nu<-3
  # theta_u<-c(rep(0,nu),c(1,2),0.5,0.1,c(3,3,3))
  grad<-function(theta_u,wpluswt,wtw,w,x,I,yo,nu,no,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p means the dimension of X.
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)] #transform of sigma2
    lamda<-theta_u[(nu+p+2)] #transform of sigma rho
    
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    
    
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
    
    #############################- M Direct inverse
    # tic()
    # dlogdetM_by_dlamda<-sum(diag(solve(M)%*%dM_by_dlamda))
    # print(dlogdetM_by_dlamda)
    # toc()
    
    # print(class(M))
    # print(class(dM_by_dlamda))
    #############################- M Choleskey inverse
    # tic()
    dlogdetM_by_dlamda<-sum(diag(solve(Cholesky(M),dM_by_dlamda)))
    # print(dlogdetM_by_dlamda)
    # toc()
    ########################
    
    div_h_lamda<-0.5*dlogdetM_by_dlamda-0.5*exp(-gama)*t(y-x%*%beta)%*%dM_by_dlamda%*%(y-x%*%beta)-lamda/omega_lamda
    
    
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_gamma),as.vector(div_h_lamda))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  # nu<-3
  # theta_u<-c(rep(1,3),-1,-2,0.1,0.2,rep(9,3))
  
  cal_Lbound_Lo<-function(theta,no,mu,d,B,xo,yo,wpluswt,wtw){
    
    p<-ncol(x) #length of betas
    
    beta<-theta[(1):(p)]
    gama<-theta[(p+1)]
    lamda<-theta[(p+2)]
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    
    p_<-ncol(B)
    n<-length(y)
    M<-I-rho*wpluswt+rho^2*wtw
    
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_psi=10000
    omega_lamda=10000
    
    D<-Diagonal(n=length(d),x=as.vector(d)) #making D sparse and diagonal
    
    
    
    M=I-rho*wpluswt+rho*rho*wtw
    M_oo=M[1:no,1:no]
    M_uu=M[(no+1):n,(no+1):n]
    M_ou=M[1:no,(no+1):n]
    M_uo=t(M_ou)
    
    V_oo_inv=M_oo-M_ou%*%solve(M_uu,M_uo)
    V_oo_inv=forceSymmetric(V_oo_inv)
    
    log_det<-determinant(V_oo_inv)[[1]][[1]]
    loglike2<--0.5*no*log(2*pi)-0.5*no*gama+
      0.5*log_det-0.5*exp(-gama)*t(yo-xo%*%beta)%*%V_oo_inv%*%(yo-xo%*%beta) # p_y efficient
    
    loglike3<-dmvnorm(x=as.vector(beta),mean = rep(0,ncol(x)),sigma=omega_beta*diag(ncol(x)),
                      log=T,checkSymmetry = F) #p_beta
    loglike4<-dnorm(x=gama,mean = 0,sd=omega_gamma) #p_gama (variance)
    
    loglike5<-dnorm(x=lamda,mean = 0,sd=omega_lamda) #p_lamda (variance)
    
    # logq<-dmvnorm(x=as.vector(theta_u),mean = as.vector(mu),sigma =B%*%t(B)+D^2,log = T,checkSymmetry = F) # q1
    
    ################################################################################ Q using Davids method
    
    r<-theta-mu
    B_tilda<-Diagonal(x = (1/d))%*%B
    r_tilde<- Diagonal(x = (1/d))%*%r
    
    epsilon_hat<-t(B_tilda)%*%r_tilde
    I_BtB<-forceSymmetric(Diagonal(n=p_,x=1)+t(B_tilda)%*%B_tilda) # o/w, i.e if not "forceSymmetric" logq may generate NaN values 
    
    ######## chol
    R<-as(chol(I_BtB),"CsparseMatrix")
    r_tilde2<-solve(t(R),epsilon_hat)
    rtcovr<-t(r_tilde)%*%r_tilde-t(r_tilde2)%*%r_tilde2 #Davids
    logdet<-sum(log(as.vector(d^2)))+sum(log(diag(R^2)))
    
    ####### choleskey
    # R<-(Matrix::Cholesky(I_BtB)) We casnt use Cholesky since I_BtB is not sparse.
    
    
    logq<--0.5*length(theta)*log(2*pi)-0.5*logdet-0.5*rtcovr #q2
    
    
    ##########################################################################
    
    # logdet<-2*sum(log(diag(chol(B%*%t(B)+D^2))))
    # logq<--0.5*nu_theta*log(2*pi)-0.5*logdet-0.5*t(r)%*%solve(B%*%t(B)+D^2,r) #q3
    
    
    
    
    return(as.numeric(loglike2)+
             as.numeric(loglike3)+as.numeric(loglike4)+
             as.numeric(loglike5)-as.numeric(logq))
    
    
  } #to calculate lower bound
  
  # 1 calculate det(M) usigng chol(), and 2 uses choleskey(). 2 is very fast
  
  #initial values for lamdha
  #p=2
  # mu<-rep(0.1,nu)
  mu<-c(start_theta)
  d<-rep(0.1,totpara)
  B_vec<-matrix(rep(0.1,(totpara)*p),ncol=1)
  B<-matrix(as.vector(B_vec),ncol=p)
  B[upper.tri(B)] <- 0
  
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
  
  ###################################### can be removed
  RMS_delta_mu<-rep(0,totpara)
  RMS_delta_B<-rep(0,totpara*p)
  RMS_delta_d<-rep(0,totpara)
  ####################################
  
  # creating a vector for storing lowerbounds
  
  Lbound_Lo<-c()
  alllamda<-c() #To store lamda in each iteration
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) 
  
  
  #i=1
  for(i in 1:N){
    
    zeta<-(rnorm(n=(p+totpara)))
    z<-matrix(zeta[1:p],ncol = 1)
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    
    
    # generate thetas
    theta<-mu+B%*%z+d*as.vector(epsilon)
    D<-diag(d)
    # generate yus from p(yu | yo, theta)
    gamma.lamda<-theta[-(1:rc)]
    rho<-(exp(gamma.lamda[2])-1)/(exp(gamma.lamda[2])+1)
    sigma2y<-exp(gamma.lamda[1])
    beta<-theta[(1:rc)]
    
    M=I-rho*wpluswt+rho*rho*wtw
    
    
    M_oo=M[1:no,1:no]
    M_uu=M[(no+1):n,(no+1):n]
    M_ou=M[1:no,(no+1):n]
    M_uo=t(M_ou)
    
    
    cholMuu<-Cholesky(M_uu)
    mu_ugo<-xu%*%beta-solve(cholMuu,M_uo)%*%(yo-xo%*%beta) # efficient cal of con mean
    # mu_ugo<-xu%*%beta+V_uo%*%solve(V_oo,yo-xo%*%beta)
    sigma_ugo<-sigma2y*solve(cholMuu,I_nu) # efficient cal of con var
    # sigma_ugo<-sigma2y*(V_uu-V_uo%*%solve(V_oo,V_ou))
    
    # yu.star=rmvnorm(1,as.vector(mu_ugo),as.matrix(sigma_ugo), #simulate from mvtnorm pkg
    #                 checkSymmetry = F) 
    yu.star<-rmvn(1,as.vector(mu_ugo),as.matrix(sigma_ugo)) #simulate from mvnfast pkg, faster
    
    theta_u<-c(as.vector(yu.star),as.vector(theta))
    
    gradg_theta<-grad(theta_u=theta_u,w=w,x=x,I=I,yo=yo,
                      wpluswt=wpluswt,wtw=wtw,nu=nu,no=no,n=n) # gratient of log h(theta)
    
    
    # ###-part-2 cal1
    # tic()
    # part2<-solve(forceSymmetric((B%*%t(B)+D^2)))%*%(B%*%z+d*as.vector(epsilon))# the second part inside the expectation
    # toc()                # Direct calculation
    # print(part2[2:4])
    
    ###-part-2 cal2
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    # 
    #Construct UB estimates
    grad_mu<-gradg_theta-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradg_theta))
    grad_B<-t(t(z)%x%I_theta)%*%(gradg_theta-dlogq_by_dtheta)
    gradd<-epsilon*(gradg_theta-dlogq_by_dtheta)
    grad_d<-matrix(gradd,ncol = 1)
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    # RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    # RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    # RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B_vec<-B_vec+rho_B*grad_B
    B<-matrix(as.vector(B_vec),ncol=p)
    B[upper.tri(B)] <- 0
    
    d<-d+rho_d*grad_d
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    ######################################## can be removed
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    #############################################
    
    # calculate lowerbound
    
    # # tic()
    Lbound_Lo[i]<-0
    # # toc()
    all_paras[i,]<-as.vector(mu)
    # Lbound[i]<-0
    ##############################
  }
  
  return(list(mu=mu,B=B,d=d,Lbound_Lo=Lbound_Lo,
              all_paras=all_paras,predictyu.draw=yu.star,
              post.mean=mu_ugo,post.cov=sigma_ugo))
} # HVB-NoB



SEM_MAR_VB_aug_ADADELTA_Blocking<-function(x,N,N1,w,yo,no,
                                           p,start_theta,bsize,yu.start){ 
  # p=1
  x<-cbind(1,x)
  rc=ncol(x) #number of regression covariates. i.e. betas
  lc=0     #number of logistic covariates. i.e. betas
  totpara<-rc+lc+2
  adapt_epsilon=10^(-6)
  v=0.95
  I_p<-Diagonal(p)
  
  n<-nrow(w)
  # nu<-sum(is.na(y_ou))
  nu<-n-no
  I<-Diagonal(n)
  I_nu<-Diagonal(nu)
  missingindexes<-(no+1):n
  xu<-x[missingindexes,]
  xo<-x[setdiff(1:n,missingindexes),]
  
  
  # bsize=10
  k=floor(nu/bsize) # # of full blocks
  reminder<-nu%%bsize # length of remaining incomplte block
  
  print(paste("# of full blocks:", k))
  print(paste("length of incomplte block:", reminder))
  
  wpluswt=(w+t(w))
  wtw=t(w)%*%(w)
  
  # nu<-3
  # theta_u<-c(rep(0,nu),c(1,2),0.5,0.1,c(3,3,3))
  grad<-function(theta_u,wpluswt,wtw,w,x,I,yo,nu,no,n){
    # m<-ncol(x)
    p<-ncol(x) # Here p means the dimension of X.
    yu<-theta_u[1:nu]
    beta<-theta_u[(nu+1):(nu+p)]
    gama<-theta_u[(nu+p+1)] #transform of sigma2
    lamda<-theta_u[(nu+p+2)] #transform of sigma rho
    
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_lamda=10000
    
    
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
    
    #############################- M Direct inverse
    # tic()
    # dlogdetM_by_dlamda<-sum(diag(solve(M)%*%dM_by_dlamda))
    # print(dlogdetM_by_dlamda)
    # toc()
    
    # print(class(M))
    # print(class(dM_by_dlamda))
    #############################- M Choleskey inverse
    # tic()
    dlogdetM_by_dlamda<-sum(diag(solve(Cholesky(M),dM_by_dlamda)))
    # print(dlogdetM_by_dlamda)
    # toc()
    ########################
    
    div_h_lamda<-0.5*dlogdetM_by_dlamda-0.5*exp(-gama)*t(y-x%*%beta)%*%dM_by_dlamda%*%(y-x%*%beta)-lamda/omega_lamda
    
    
    
    fullgrad<-c(as.vector(div_h_beta),
                as.vector(div_h_gamma),as.vector(div_h_lamda))
    
    return(fullgrad)
    
  } #Calculate gradient of original parameters, theta
  
  # nu<-3
  # theta_u<-c(rep(1,3),-1,-2,0.1,0.2,rep(9,3))
  cal_Lbound_Lo<-function(theta,no,mu,d,B,xo,yo,
                          wpluswt,wtw){
    
    p<-ncol(x) #length of betas
    
    beta<-theta[(1):(p)]
    gama<-theta[(p+1)]
    lamda<-theta[(p+2)]
    rho<-(exp(lamda)-1)/(exp(lamda)+1)
    
    
    p_<-ncol(B)
    n<-length(y)
    M<-I-rho*wpluswt+rho^2*wtw
    
    
    #Hyper parameters
    omega_beta=10000 # prior of betas is normal with zero mean and large variance
    omega_gamma=10000
    omega_psi=10000
    omega_lamda=10000
    
    D<-Diagonal(n=length(d),x=as.vector(d)) #making D sparse and diagonal
    
    
    
    M=I-rho*wpluswt+rho*rho*wtw
    M_oo=M[1:no,1:no]
    M_uu=M[(no+1):n,(no+1):n]
    M_ou=M[1:no,(no+1):n]
    M_uo=t(M_ou)
    
    V_oo_inv=M_oo-M_ou%*%solve(M_uu,M_uo)
    V_oo_inv=forceSymmetric(V_oo_inv)
    
    log_det<-determinant(V_oo_inv)[[1]][[1]]
    loglike2<--0.5*no*log(2*pi)-0.5*no*gama+
      0.5*log_det-0.5*exp(-gama)*t(yo-xo%*%beta)%*%V_oo_inv%*%(yo-xo%*%beta) # p_y efficient
    
    loglike3<-dmvnorm(x=as.vector(beta),mean = rep(0,ncol(x)),sigma=omega_beta*diag(ncol(x)),
                      log=T,checkSymmetry = F) #p_beta
    loglike4<-dnorm(x=gama,mean = 0,sd=omega_gamma) #p_gama (variance)
    
    loglike5<-dnorm(x=lamda,mean = 0,sd=omega_lamda) #p_lamda (variance)
    
    # logq<-dmvnorm(x=as.vector(theta_u),mean = as.vector(mu),sigma =B%*%t(B)+D^2,log = T,checkSymmetry = F) # q1
    
    ################################################################################ Q using Davids method
    
    r<-theta-mu
    B_tilda<-Diagonal(x = (1/d))%*%B
    r_tilde<- Diagonal(x = (1/d))%*%r
    
    epsilon_hat<-t(B_tilda)%*%r_tilde
    I_BtB<-forceSymmetric(Diagonal(n=p_,x=1)+t(B_tilda)%*%B_tilda) # o/w, i.e if not "forceSymmetric" logq may generate NaN values 
    
    ######## chol
    R<-as(chol(I_BtB),"CsparseMatrix")
    r_tilde2<-solve(t(R),epsilon_hat)
    rtcovr<-t(r_tilde)%*%r_tilde-t(r_tilde2)%*%r_tilde2 #Davids
    logdet<-sum(log(as.vector(d^2)))+sum(log(diag(R^2)))
    
    ####### choleskey
    # R<-(Matrix::Cholesky(I_BtB)) We casnt use Cholesky since I_BtB is not sparse.
    
    
    logq<--0.5*length(theta)*log(2*pi)-0.5*logdet-0.5*rtcovr #q2
    
    
    ##########################################################################
    
    # logdet<-2*sum(log(diag(chol(B%*%t(B)+D^2))))
    # logq<--0.5*nu_theta*log(2*pi)-0.5*logdet-0.5*t(r)%*%solve(B%*%t(B)+D^2,r) #q3
    
    
    
    
    return(as.numeric(loglike2)+
             as.numeric(loglike3)+as.numeric(loglike4)+
             as.numeric(loglike5)-as.numeric(logq))
    
    
  } #to calculate lower bound
  
  # 1 calculate det(M) usigng chol(), and 2 uses choleskey(). 2 is very fast
  make.w_x_<-function(w,x,k,reminder){
    
    if(reminder==0){
      
      w.list<-list()
      x.list<-list()
      
      
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
        
        
        
        # m<-rev(m)
        current<-current+bsize
        # print("no half block")
      }}else{
        
        w.list<-list()
        x.list<-list()
        
        
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
            
            
          }
        }
        # print("There is a half block")
      }
    
    
    return(list(w.list=w.list,x.list=x.list))
  }
  
  sim_yugyotheta_<-function(yo,listw_x_m,    # After adding un-balanced # of blocks 
                            yu.start,beta,
                            sigma2y,rho,
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
          
          current_state <- proposal_state
          
          
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
          current_state <- proposal_state
          accepts<-accepts+1
          
          
        }
      }
      
      
      
      
      yu.chain[i,] <- current_state
    }
    # toc()
    
    
    
    return(list(yugyom.chain=yu.chain,accepts=accepts,
                co.dis=co.dis))
    
  }
  #initial values for lamdha
  #p=2
  # mu<-rep(0.1,nu)
  mu<-c(start_theta)
  d<-rep(0.1,totpara)
  B_vec<-matrix(rep(0.1,(totpara)*p),ncol=1)
  B<-matrix(as.vector(B_vec),ncol=p)
  B[upper.tri(B)] <- 0
  
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
  
  ###################################### can be removed
  RMS_delta_mu<-rep(0,totpara)
  RMS_delta_B<-rep(0,totpara*p)
  RMS_delta_d<-rep(0,totpara)
  ####################################
  
  # creating a vector for storing lowerbounds
  
  Lbound_Lo<-c()
  alllamda<-c() #To store lamda in each iteration
  all_paras<-matrix(rep(0,N*length(mu)),nrow = N) 
  
  listw_x_m<-make.w_x_(w=w,x=x,k=k,reminder =reminder)
  #i=1
  for(i in 1:N){
    
    zeta<-(rnorm(n=(p+totpara)))
    z<-matrix(zeta[1:p],ncol = 1)
    epsilon<-matrix(zeta[p+1:(length(zeta)-p)],ncol=1)
    d<-as.vector(d) # d upated from algorithm is a matrix, we need to change to a vector
    
    
    # generate thetas
    theta<-mu+B%*%z+d*as.vector(epsilon)
    D<-diag(d)
    # generate yus from p(yu | yo, theta)
    gamma.lamda<-theta[-(1:rc)]
    rho<-(exp(gamma.lamda[2])-1)/(exp(gamma.lamda[2])+1)
    sigma2y<-exp(gamma.lamda[1])
    beta<-theta[(1:rc)]
    
    
    
    # tic()
    # M=I-rho*wpluswt+rho*rho*wtw
    # M_oo=M[1:no,1:no]
    # M_uu=M[(no+1):n,(no+1):n]
    # M_ou=M[1:no,(no+1):n]
    # M_uo=t(M_ou)
    # cholMuu<-Cholesky(M_uu)
    # mu_ugo<-xu%*%beta-solve(cholMuu,M_uo)%*%(yo-xo%*%beta) # efficient cal of con mean
    # sigma_ugo<-sigma2y*solve(cholMuu,I_nu)
    # yu.star<-rmvn(1,as.vector(mu_ugo),as.matrix(sigma_ugo)) #simulate from mvnfast pkg, faster
    # toc()
    
    # tic()
    yu.star.result<-sim_yugyotheta_(yo=yo,listw_x_m=listw_x_m,
                                    yu.start=yu.start,beta=beta,
                                    sigma2y=sigma2y,rho=rho
                                    ,N1=N1,k=k)
    yu.star=yu.star.result$yugyom.chain[N1,]
    yu.start=yu.star # from the second VB iteration, the starting value of MCMC
    # is set to previous yu
    # toc()
    
    theta_u<-c(as.vector(yu.star),as.vector(theta))
    
    gradg_theta<-grad(theta_u=theta_u,w=w,x=x,I=I,yo=yo,
                      wpluswt=wpluswt,wtw=wtw,nu=nu,no=no,n=n) # gratient of log h(theta)
    
    
    # ###-part-2 cal1
    # tic()
    # part2<-solve(forceSymmetric((B%*%t(B)+D^2)))%*%(B%*%z+d*as.vector(epsilon))# the second part inside the expectation
    # toc()                # Direct calculation
    # print(part2[2:4])
    
    ###-part-2 cal2
    
    D_t2_inv<-Diagonal(n=length(d),x=(1/d^2)) # calculate inverse using woodbury formula
    part2<-(D_t2_inv-D_t2_inv%*%B%*%solve(I_p+t(B)%*%D_t2_inv%*%B)%*%t(B)%*%D_t2_inv)%*%(B%*%z+d*as.vector(epsilon))
    dlogq_by_dtheta<--part2
    # 
    #Construct UB estimates
    grad_mu<-gradg_theta-dlogq_by_dtheta
    I_theta<-Diagonal(n=length(gradg_theta))
    grad_B<-t(t(z)%x%I_theta)%*%(gradg_theta-dlogq_by_dtheta)
    gradd<-epsilon*(gradg_theta-dlogq_by_dtheta)
    grad_d<-matrix(gradd,ncol = 1)
    
    # calculate new learning rates#################
    E_g2_mu<-v*E_g2_mu+(1-v)*grad_mu^2
    E_g2_B<-v*E_g2_B+(1-v)*grad_B^2
    E_g2_d<-v*E_g2_d+(1-v)*grad_d^2
    
    RMS_g_mu<-sqrt(E_g2_mu+adapt_epsilon)
    RMS_g_d<-sqrt(E_g2_d+adapt_epsilon)
    RMS_g_B<-sqrt(E_g2_B+adapt_epsilon)
    
    # RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    # RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    # RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    
    
    rho_mu<-(RMS_delta_mu/RMS_g_mu)
    rho_B<-(RMS_delta_B/RMS_g_B)
    rho_d<-(RMS_delta_d/RMS_g_d)
    
    mu<-mu+rho_mu*grad_mu
    B_vec<-B_vec+rho_B*grad_B
    B<-matrix(as.vector(B_vec),ncol=p)
    B[upper.tri(B)] <- 0
    
    d<-d+rho_d*grad_d
    
    # calculate E_delta2_mu for this iteration to use in next iteration
    E_delta2_mu<-v*E_delta2_mu+(1-v)*(rho_mu*grad_mu)^2
    E_delta2_B<-v*E_delta2_B+(1-v)*(rho_B*grad_B)^2
    E_delta2_d<-v*E_delta2_d+(1-v)*(rho_d*grad_d)^2
    
    ######################################## can be removed
    RMS_delta_mu<-sqrt(E_delta2_mu+adapt_epsilon)
    RMS_delta_B<-sqrt(E_delta2_B+adapt_epsilon)
    RMS_delta_d<-sqrt(E_delta2_d+adapt_epsilon)
    #############################################
    
    # calculate lowerbound
    
    # # tic()
    Lbound_Lo[i]<-0
    # # toc()
    all_paras[i,]<-as.vector(mu)
    # Lbound[i]<-0
    ##############################
  }
  
  return(list(mu=mu,B=B,d=d,Lbound_Lo=Lbound_Lo,
              all_paras=all_paras,predictyu.draw=yu.star))
} ## HVB-G

Marginallikemethod<-function(x,y,w,no,model=1){
  
  n=ncol(w)
  x=cbind(rep(1,n),x)
  nu<-n-no
  I=Diagonal(n)
  p=ncol(x)
  
  
  # class(model)
  wpluswt=(w+t(w))
  wwt=t(w)%*%(w)
  
  marginlikelihood<-function(rho,x,y,w,wpluswt=wpluswt,wwt=wwt,no=no,I=I,model=model){
    
    n=nrow(w)
    
    
    ###calculate vss
    AtA=(I-rho*(wpluswt)+rho*rho*wwt)
    v=solve(AtA)
    v=forceSymmetric(v) # make symetric 
    voo=v[1:no,1:no]
    
    A=I-rho*w
    
    if(model==1){
      x_till=x
    }else{
      x_till=solve(A,x)
    }
    
    x_tillo=x_till[(1:no),]
    yo=y[1:no]
    
    ### beta_hat and omega_hat
    
    betahat=solve(t(x_tillo)%*%solve(voo,x_tillo))%*%t(x_tillo)%*%solve(voo,yo)
    
    #3
    ro=(yo-x_tillo%*%betahat)
    #4
    omegahat=(t(ro)%*%solve(voo,ro))/no
    
    log_det=-sum(log(diag(chol(voo))))
    
    my_loglik=(-no/2)*log(2*pi)-(no/2)*log(omegahat)+log_det-(no/2)
    return(-as.numeric(my_loglik))
    
  }
  
  # rho0=0.1
  res<-optimize(f=marginlikelihood,lower=-0.9999,upper=0.9999,
                x=x,y=y,w=w,wpluswt=wpluswt,wwt=wwt,no=no,I=I,model=model)
  rho<-res$minimum
  
  AtA=(I-rho*(wpluswt)+rho*rho*wwt)
  v=solve(AtA)
  v=forceSymmetric(v) 
  voo=v[1:no,1:no]
  
  A=I-rho*w
  
  if(model==1){
    x_till=x
  }else{
    x_till=solve(A,x)
  }
  
  x_tillo=x_till[(1:no),]
  yo=y[(1:no)]
  
  
  betahat=solve(t(x_tillo)%*%solve(voo,x_tillo))%*%t(x_tillo)%*%solve(voo,yo)
  ro=(yo-x_tillo%*%betahat)
  omegahat=(t(ro)%*%solve(voo,ro))/no
  
  l=list("rho"=rho,"beta"=betahat,"omegao"=omegahat)
  return(l)
  
  
} # MML- (Suesse, 2018)


