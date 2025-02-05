source("source.R")




######################### data



######################### data



sqrtn<-25
g<-graph.lattice(c(sqrtn,sqrtn)) 
g<-connect.neighborhood(g,1) #connect all vertices by one edge.
wmat<-as_adjacency_matrix(g)
w<-mat2listw(wmat,style="W")
w=listw2mat(w)
w<-as(w,"CsparseMatrix")


p=10


para=sample(1:5,(p+1),replace = T)
para=c(para,c(0.8,1))

sim=simulateSEM(para=para,p=p,weightmat=w)
xsim<-sim$Independent
ysim<-sim$Dependent

x=as.matrix(xsim)
y=ysim

p_sample<-0.25
splitted_sample=splitter(x=x,y=y,w=w,p_sample=p_sample)


x.split=splitted_sample$X
y.split=splitted_sample$Y
w.split=splitted_sample$W
no=splitted_sample$ns

n<-ncol(w)
nu=n-no

yu<-y.split[(no+1):n]
y_ou<-c(y.split[1:no],rep(NA,nu)) # for HMC
yo<-y.split[1:no]
############## starting values



n<-ncol(w)
no<-n-nu

# For thetas
reg.model<-estimate_regression(x=x.split,yo=yo,no=no)
reg.model.summary<-summary(reg.model)
sigma2y0<-(reg.model.summary$sigma)^2
gama0<-log(sigma2y0)
beta0<-as.numeric((reg.model.summary$coefficients)[,1])
lamda0<-log(1+0.001)-log(1-0.001)
paras.st<-c(beta0,gama0,lamda0)

# for missing values, yu


yu0<-sim_yugyo(yo,x.split,w.split,paras.st)
plot(yu,yu0)
paras.st<-c(yu0,paras.st)




#################################### Run Algorithms

#JVB

N=10000

tic()
fit_SEM_MAR_VB_thetayu<-SEM_MAR_VB.new(x=x.split,N=N,
                                       w=w.split,no=no,y=y.split,
                                       p=4,start_theta=paras.st)
toc()

# HVB-NoB


n<-ncol(w)
no<-n-nu
reg.model<-estimate_regression(x=x.split,yo=yo,no=no)
reg.model.summary<-summary(reg.model)
sigma2y0<-(reg.model.summary$sigma)^2
gama0<-log(sigma2y0)
beta0<-as.numeric((reg.model.summary$coefficients)[,1])
lamda0<-log(1+0.001)-log(1-0.001)
paras.st<-c(beta0,gama0,lamda0)
paras.st<-c(paras.st)


yo<-y.split[1:no]
tic()
fit_SEM_MAR_VB_thetaAug<-SEM_MAR_VB_aug_ADADELTA(x=x.split,N=N,
                                                 w=w.split,no=no,yo=yo,p=4,
                                                 start_theta=paras.st)
toc()




# HMC


iter=10000

code.path='SEM_MAR_HMC.stan'

data.list=make.SATN.data(x=x.split,N=N,w=w.split,y=y_ou,no=no)


tic()
fit_SEM_HMC_thetayu<-stan(file = code.path,
                          data = data.list,chains = 1,iter = iter) # default warm upsize=iter/2.
toc()

#Extract MCMC chains from HMC output
mcmc.chains.T<-as.matrix(fit_SEM_HMC_thetayu)
HMC.theta.chains<-mcmc.chains.T[,(1:(ncol(mcmc.chains.T)-nu-1))]
yu.chains<-mcmc.chains.T[,-(1:(ncol(mcmc.chains.T)-nu-1))]
HMC.yu.chains<-yu.chains[,-(nu+1)]


# Marginal MLE

MML.estimates<-Marginallikemethod(x.split,y.split,w.split,no)


###### Sample from posteriors of theta

# JVB
sample.size<-nrow(HMC.theta.chains)
VB1.sample<-matrix(rep(0,sample.size*(13+nu)),ncol = 13+nu)

for(i in 1:sample.size){
  VB1.sample[i,]<-simulate.VB.MAR_vb1(mu=fit_SEM_MAR_VB_thetayu$mu,
                                      B=fit_SEM_MAR_VB_thetayu$B,
                                      d=fit_SEM_MAR_VB_thetayu$d,nu=nu,x=x.split,p=4)
  
}

round(apply(VB1.sample,2,mean)[-(1:nu)],4)
round(apply(VB1.sample,2,sd)[-(1:nu)],4)

### HVB


VB2.sample<-matrix(rep(0,sample.size*(13)),ncol = 13)

for(i in 1:sample.size){
  VB2.sample[i,]<-simulate.VB.SEM_vb2(mu=fit_SEM_MAR_VB_thetaAug$mu,
                                      B=fit_SEM_MAR_VB_thetaAug$B,
                                      d=fit_SEM_MAR_VB_thetaAug$d,x=x,p=4)
  
}
round(apply(VB2.sample,2,mean),4)
round(apply(VB2.sample,2,sd),4)


############# Density plots
# B0





name.list<-c("HMC","VB-yutheta","VB-upall")

values=c("HMC"= "brown4","VB-yutheta"="blue","VB-upall"="green")


# beta0

df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(HMC.theta.chains[,"beta[1]"], VB1.sample[,(nu+1)],VB2.sample[,1])
)


pp1 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 5)) +
  scale_y_continuous(limits = c(0, 4)) +
  scale_color_manual(values = values, breaks = c("HMC", "VB-yutheta","VB-upall"),
                     labels = c("HMC", "JVB", "HVB-NoB")) +  # Specify breaks
  labs(
    x = expression(beta[0]),
    y = "Density",
    color = "") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.27, 0.86),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 22),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )+geom_vline(xintercept = para[1], linetype = "dashed", color = "orange")

pp1

# beta5


df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(HMC.theta.chains[,"beta[6]"], VB1.sample[,(nu+6)],VB2.sample[,6])
)


ppb5 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 5)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  scale_color_manual(values = values, breaks = c("HMC", "VB-yutheta","VB-upall"),
                     labels = c("HMC", "VB-1", "VB-2")) +  # Specify breaks
  labs(
    x = expression(beta[6]),
    y = "Density",
    color = "") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.86),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 22),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )+geom_vline(xintercept = para[6], linetype = "dashed", color = "orange")

ppb5

# beta10


df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(HMC.theta.chains[,"beta[11]"], VB1.sample[,(nu+11)],VB2.sample[,11])
)


ppb10 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 5)) +
  scale_y_continuous(limits = c(0, 6.75)) +
  scale_color_manual(values = values, breaks = c("HMC", "VB-yutheta","VB-upall"),
                     labels = c("HMC", "VB-1", "VB-2")) +  # Specify breaks
  labs(
    x = expression(beta[10]),
    y = "Density",
    color = "") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.86),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 22),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )+geom_vline(xintercept = para[11], linetype = "dashed", color = "orange")

ppb10

# beta8


df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(HMC.theta.chains[,"beta[9]"], VB1.sample[,(nu+9)],VB2.sample[,9])
)


ppb8 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(1, 5)) +
  scale_y_continuous(limits = c(0, 4.5)) +
  scale_color_manual(values = values, breaks = c("HMC", "VB-yutheta","VB-upall"),
                     labels = c("HMC", "VB-1", "VB-2")) +  # Specify breaks
  labs(
    x = expression(beta[8]),
    y = "Density",
    color = "") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.86),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 22),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )+geom_vline(xintercept = para[9], linetype = "dashed", color = "orange")

ppb8

#signa2y
para
hmc.sigma2y.smple<-exp(HMC.theta.chains[,"gama"])
df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(hmc.sigma2y.smple, VB1.sample[,(nu+p+2)],VB2.sample[,p+2])
)

pp2 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 3)) +
  scale_y_continuous(limits = c(0, 3.5)) +
  scale_color_manual(values = values, breaks = c("HMC", "VB-yutheta","VB-upall"),
                     labels = c("HMC", "VB-1", "VB-2")) +  # Specify breaks
  labs(
    x = expression(sigma[e]^2),
    y = "Density",
    color = "") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.25, 0.72),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 22),  # Adjust the size as needed
    legend.text = element_text(size = 5))+geom_vline(xintercept = para[13], 
                                                     linetype = "dashed", color = "orange")

pp2

# rho

## rho

hmc.rho.smple<-((exp(HMC.theta.chains[,"lamda"]) -1)/(exp(HMC.theta.chains[,"lamda"]) +1))
df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(hmc.rho.smple, VB1.sample[,(nu+p+3)],VB2.sample[,p+3])
)

pp3 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(-0.1, 1)) +
  scale_y_continuous(limits = c(0, 8.6)) +
  scale_color_manual(values = values, breaks = c("HMC", "VB-yutheta","VB-upall")) +  # Specify breaks
  labs(
    x = expression(rho),
    y = "Density",
    color = "") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.25, 0.72),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 22),  # Adjust the size as needed
    legend.text = element_text(size = 5))+geom_vline(xintercept = para[12], 
                                                     linetype = "dashed", color = "orange")

pp3


pp=pp1+ppb5+ppb8+ppb10+pp2+pp3+plot_layout(ncol = 3,nrow = 2)

########## missing values

####### Posterior mean and sd comparison of Missing values

tic()
HVB.post.sample.yu<-apply(VB2.sample, 1, post.sample.generate,yo=yo,x=x.split,w=w.split)
toc()


# VB Vs HMC means

data <- data.frame(
  x = apply(HMC.yu.chains,2,mean),  # Replace with your x values
  vb1 = apply(VB1.sample[,1:nu],2,mean),  # Replace with your y values
  vb2=apply(HVB.post.sample.yu, 1, mean) # we have posterior mean
)

pp1 <- ggplot(data, aes(x = x)) +
  geom_point(aes(y = vb1, color = "vb1"), size = 2) +
  geom_point(aes(y = vb2, color = "VB2"), size = 2) +  # Add vb2 values
  labs(
    title = "",
    x = "HMC posterior means",
    y = "VB posterior means"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c("vb1" = "blue", "VB2" = "green"),  # Adjust colors as needed
    labels = c("JVB", "HVB-NoB")  # Adjust labels as needed
  ) +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = c(0.25,0.9),  # Position the legend at the top of the plot
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a border around the plot
    
    axis.text = element_text(size = 22),  # Adjust the size as needed
    axis.title = element_text(size = 22),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5, size = 25),  # Adjust the size as needed
    legend.text = element_text(size = 22)  # Adjust the size as needed
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "green") +
  guides(color = guide_legend(title = ""))  # legend title

pp1  # Display the plot


# VB Vs HMC sd


data <- data.frame(
  x = apply(HMC.yu.chains,2,sd),  # Replace with your x values
  vb1 = apply(VB1.sample[,1:nu],2,sd),  # Replace with your y values
  vb2=apply(HVB.post.sample.yu, 1, sd)
)
pp2 <- ggplot(data, aes(x = x)) +
  geom_point(aes(y = vb1, color = "vb1"), size = 2) +
  geom_point(aes(y = vb2, color = "VB2"), size = 2) +  # Add vb2 values
  labs(
    title = "",
    x = "HMC posterior standard deviations",
    y = "VB posterior standard deviations"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c("vb1" = "blue", "VB2" = "green"),  # Adjust colors as needed
    labels = c("JVB", "HVB-NoB")  # Adjust labels as needed
  ) +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = c(0.8,0.25),  # Position the legend at the top of the plot
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a border around the plot
    
    axis.text = element_text(size = 22),  # Adjust the size as needed
    axis.title = element_text(size = 22),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5, size = 25),  # Adjust the size as needed
    legend.text = element_text(size = 22)  # Adjust the size as needed
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "green") +
  guides(color = guide_legend(title = ""))  # legend title

pp2  # Display the plot


##### VB1 JVB-lower bound


# Convert time series to dataframe
ts_data <- data.frame(Iteration = time(fit_SEM_MAR_VB_thetayu$Lbound), 
                      Value = fit_SEM_MAR_VB_thetayu$Lbound)

# Plot using ggplot2
pp.con.vb1=ggplot(ts_data, aes(x = Iteration, y = Value)) +
  geom_line(color = "blue",size = 0.25) +
  labs(x = "Iteration", y = "Lower bound of JVB")+
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.29, 0.82),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )+ggtitle("JVB")

pp.con.vb1
####VB-2

# NoB

vb2.trojectry=(fit_SEM_MAR_VB_thetaAug$all_paras)[,c(1,2,12,13)]
vb2.trojectry[,3]<-exp(vb2.trojectry[,3])
vb2.trojectry[,4]<-((exp(vb2.trojectry[,4]) -1)/(exp(vb2.trojectry[,4]) +1))
# vb1.trojectry[10000,]


vb2.trojectry<-as.data.frame(vb2.trojectry)
vb2.trojectry$Time<-1:10000
vb1.trojectry_melted <- melt(vb2.trojectry, id.vars = "Time")

legend_labels <- c(expression(beta[0]),expression(beta[1]),
                   expression(sigma[e]^2), expression(rho))

# Define line types
line_types <- c("solid", "dashed", "dotted", "dotdash")

pp.con.vbnob <- ggplot(data = vb1.trojectry_melted, aes(x = Time, y = value, color = variable, linetype = variable)) +
  geom_line(size = 0.25) +
  scale_color_manual(values = rep("green", length(legend_labels)), name = NULL, labels = legend_labels) +
  scale_linetype_manual(values = line_types, name = NULL, labels = legend_labels) + # Specify line types and remove linetype legend title
  labs(title = "",
       x = "Iteration", y = "Value") +
  theme_minimal() +
  theme(
    
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.85),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )+ggtitle("HVB-NoB")+ylim(0,4.5)

pp.con.vbnob


pp=pp.con.vb1+pp.con.vbnob+plot_layout(ncol = 2,nrow = 1)

## HMC trace plots

beta0<-mcmc.chains.T[,"beta[1]"]
beta1<-mcmc.chains.T[,"beta[2]"]
gama<-mcmc.chains.T[,"gama"]
sigma2y<-exp(gama)
lamda<-mcmc.chains.T[,"lamda"]
rho<-((exp(lamda) -1)/(exp(lamda) +1))

###


beta0_df <- data.frame(
  Time = time(beta0),
  Value = as.numeric(beta0)
)

hmc.betao=ggplot(beta0_df, aes(x = Time, y = Value)) +
  geom_line(color="brown4") +
  labs(title = expression(beta[0]),
       x = "Iteration",
       y = "Value") +
  theme(
    
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.8),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )



beta1_df <- data.frame(
  Time = time(beta1),
  Value = as.numeric(beta1)
)

hmc.beta1=ggplot(beta1_df, aes(x = Time, y = Value)) +
  geom_line(color="brown4") +
  labs(title = expression(beta[1]),
       x = "Iteration",
       y = "Value") +
  theme(
    
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.8),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )


sigma2y_df <- data.frame(
  Time = time(sigma2y),
  Value = as.numeric(sigma2y)
)

hmc.sigma2y=ggplot(sigma2y_df, aes(x = Time, y = Value)) +
  geom_line(color="brown4") +
  labs(title = expression(sigma[e]^2),
       x = "Iteration",
       y = "Value") +
  theme(
    
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.8),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )





rho_df <- data.frame(
  Time = time(rho),
  Value = as.numeric(rho)
)

hmc.rho=ggplot(rho_df, aes(x = Time, y = Value)) +
  geom_line(color="brown4") +
  labs(title = expression(rho),
       x = "Iteration",
       y = "Value") +
  theme(
    
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.2, 0.8),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )






pp=hmc.betao+hmc.beta1+hmc.sigma2y+hmc.rho+plot_layout(ncol = 2,nrow = 2)



