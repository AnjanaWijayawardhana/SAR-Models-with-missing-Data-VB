source("source.R")



############################# Data set

# Extract W

hlw <- nb2listw(e80_queen,zero.policy = TRUE)
hsn <- listw2sn(hlw)
W<- as(hlw, "CsparseMatrix")

# Extract x and y


log.y<-log(elect80$pc_turnout)
X<-data.frame(log.edu=log(elect80$pc_college),
              log.house=log(elect80$pc_homeownership),
              income=(elect80$pc_income),
              log.edu.log.house=log(elect80$pc_college)*log(elect80$pc_homeownership),
              log.edu.income=log(elect80$pc_college)*(elect80$pc_income),
              log.house.income=log(elect80$pc_homeownership)*(elect80$pc_income),
              log.edu.log.house.income=log(elect80$pc_college)*log(elect80$pc_homeownership)*(elect80$pc_income)
)

standardize_columns <- function(mat) {
  # Calculate column means and standard deviations
  col_means <- colMeans(mat)
  col_sds <- apply(mat, 2, sd)
  
  # Standardize each column
  standardized_mat <- sweep(mat, 2, col_means, "-")
  standardized_mat <- sweep(standardized_mat, 2, col_sds, "/")
  
  return(standardized_mat)
}


Xst<-standardize_columns(X) # standardising covariates
X=Xst
head(X)


# Make missing values


p_sample<-0.25


splitted_sample=splitter(x=X,y=log.y,
                         w=W,p_sample=p_sample)
x.split=splitted_sample$X
y.split=splitted_sample$Y
w.split=splitted_sample$W
no=splitted_sample$ns
n=ncol(W)

# Set starting values, yu.start and theta.start

n<-ncol(W)
nu<-n-no
yo<-y.split[1:no]
yu<-y.split[-(1:no)]

#-starting values for thetas
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

########################### Implement VB algorithms ########################################

#HVB-G

N=15000

paras.st<-c(beta0,gama0,lamda0)
tic()
fit_SEM_MAR_VB_thetaAug_B<-SEM_MAR_VB_aug_ADADELTA_Blocking(x=x.split,N=N,N1=10,bsize=500,
                                                            w=w.split,no=no,yo=yo,p=4,
                                                            start_theta=paras.st,yu.start=yu0)
toc()

# JVB


paras.yu.st<-c(yu0,paras.st)


tic()
fit_SEM_MAR_VB_thetayu<-SEM_MAR_VB.new(x=x.split,N=N,
                                       w=w.split,no=no,y=y.split,
                                       p=4,start_theta=paras.yu.st)
toc()


# Marginal MLE  - (Suesse, 2018)

tic()
fit.marginal.SEM<-Marginallikemethod(x=x.split,y=y.split,w=w.split,no=no,model=1)
toc()


####### Sample from posteriors of theta


sample.size<-10000
VB1.sample<-matrix(rep(0,sample.size*(length(paras.st)+nu)),ncol = length(paras.st)+nu)

for(i in 1:sample.size){
  VB1.sample[i,]<-simulate.VB.MAR_vb1(mu=fit_SEM_MAR_VB_thetayu$mu,
                                      B=fit_SEM_MAR_VB_thetayu$B,
                                      d=fit_SEM_MAR_VB_thetayu$d,nu=nu,x=x.split,p=4)
  
}

round(apply(VB1.sample,2,mean)[-(1:nu)],4)
round(apply(VB1.sample,2,sd)[-(1:nu)],4)

#HVB-G

VB2.sample<-matrix(rep(0,sample.size*length(paras.st)),ncol = length(paras.st))

for(i in 1:sample.size){
  VB2.sample[i,]<-simulate.VB.SEM_vb2(mu=fit_SEM_MAR_VB_thetaAug_B$mu,
                                      B=fit_SEM_MAR_VB_thetaAug_B$B,
                                      d=fit_SEM_MAR_VB_thetaAug_B$d,x=x.split,p=4)
  
}
round(apply(VB2.sample,2,mean),4)
round(apply(VB2.sample,2,sd),4)

# fit.marginal.SEM$beta
############## simulate and plot Missing values

#Simulate HVB-G missing values

tic()
HVB.G.post.sample.yu<-apply(VB2.sample[1:10,], 1, generate_post_yu_MCMC_B,
                            w.split,x.split,yo=yo,500,yu0[1,],N1=5)
toc()


yu.post.mean.HVB.G<-apply(HVB.G.post.sample.yu, 1, mean) # cal postterior mean of missing values for HVB-G


## HMC vs VB mean

data <- data.frame(
  x = yu,  # Replace with your x values
  VB1 = apply(VB1.sample[,1:nu],2,mean),  # Replace with your VB1 values
  VB2 = yu.post.mean.HVB.G # Replace with your VB2 values
)

pp1 <- ggplot(data, aes(x = x)) +
  geom_point(aes(y = VB1, color = "VB-1"), size = 2) +
  geom_point(aes(y = VB2, color = "VB-2"), size = 2) +  # Add scatter plot for VB-2
  labs(
    title = "",
    x = "True missing values",
    y = "VB posterior means"
  ) +
  theme_minimal() +
  scale_color_manual(
    values = c("VB-1" = "blue", "VB-2" = "green"),  # Adjust colors as needed
    labels = c("JVB", "HVB-G")  # Adjust labels as needed
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.25, 0.9),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 22),
    plot.title = element_text(hjust = 0.5, size = 25),
    legend.text = element_text(size = 22)
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "green") +
  guides(color = guide_legend(title = ""))

pp1

# calculate MSEs

sum((yu-apply(VB1.sample[,1:nu],2,mean))^2) #JVB
sum((yu-yu.post.mean.HVB.G)^2) #HVB-G

#############  Generate density Plots

name.list<-c("VB-yutheta","VB-upall")

values=c("VB-yutheta"="blue","VB-upall"="green")


# beta0

df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(VB1.sample[,(nu+1)],VB2.sample[,1])
)


pp1 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5) +
  theme_minimal() +
  scale_x_continuous(limits = c(-.7, -0.5)) +
  scale_y_continuous(limits = c(0, 80)) +
  scale_color_manual(values = values, breaks = c("VB-yutheta","VB-upall"),
                     labels = c("JVB", "HVB-NoB")) +  # Specify breaks
  labs(
    x = "intercept",
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
  )+geom_vline(xintercept = fit.marginal.SEM$beta[1], linetype = "dashed", color = "orange")

pp1






# beta5


df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(VB1.sample[,(nu+3)],VB2.sample[,3])
)


ppb5 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, .5)) +
  scale_y_continuous(limits = c(0, 7)) +
  scale_color_manual(values = values, breaks = c("VB-yutheta","VB-upall"),
                     labels = c( "VB-1", "VB-2")) +  # Specify breaks
  labs(
    x = "log.prop.house",
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
  )+geom_vline(xintercept =fit.marginal.SEM$beta[3], linetype = "dashed", color = "orange")

ppb5

# log.prop.edu × income


df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c( VB1.sample[,(nu+6)],VB2.sample[,6])
)


ppb10 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(-.5, 1)) +
  scale_y_continuous(limits = c(0, 5)) +
  scale_color_manual(values = values, breaks = c("VB-yutheta","VB-upall"),
                     labels = c("VB-1", "VB-2")) +  # Specify breaks
  labs(
    x = "log.prop.edu × income",
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
  )+geom_vline(xintercept = fit.marginal.SEM$beta[6], linetype = "dashed", color = "orange")

ppb10

# log.prop.house × income


df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c( VB1.sample[,(nu+7)],VB2.sample[,7])
)


ppb8 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(-1, .5)) +
  scale_y_continuous(limits = c(0,6)) +
  scale_color_manual(values = values, breaks = c( "VB-yutheta","VB-upall"),
                     labels = c("VB-1", "VB-2")) +  # Specify breaks
  labs(
    x ="log.prop.edu × income",
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
  )+geom_vline(xintercept =fit.marginal.SEM$beta[7], linetype = "dashed", color = "orange")

ppb8

#signa2y

df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(VB1.sample[,(nu+ncol(x.split)+2)],VB2.sample[,ncol(x.split)+2])
)

pp2 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, .1)) +
  scale_y_continuous(limits = c(0, 3.5)) +
  scale_color_manual(values = values, breaks = c( "VB-yutheta","VB-upall"),
                     labels = c( "VB-1", "VB-2")) +  # Specify breaks
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
    legend.text = element_text(size = 5))+geom_vline(xintercept = as.numeric(fit.marginal.SEM$omegao), 
                                                     linetype = "dashed", color = "orange")

pp2

# rho

## rho


df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(VB1.sample[,(nu+ncol(x.split)+3)],VB2.sample[,ncol(x.split)+3])
)

pp3 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(-0.1, 1)) +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = values, breaks = c( "VB-yutheta","VB-upall")) +  # Specify breaks
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
    legend.text = element_text(size = 5))+geom_vline(xintercept = fit.marginal.SEM$rho, 
                                                     linetype = "dashed", color = "orange")

pp3


pp=pp1+ppb5+ppb8+ppb10+pp2+pp3+plot_layout(ncol = 3,nrow = 2)



##### Convergence plots 

## JVB lower bound


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


#### HVB-G

vb2.trojectry=(fit_SEM_MAR_VB_thetaAug_B$all_paras)[,c(1,3,9,10)]
vb2.trojectry[,3]<-exp(vb2.trojectry[,3])
vb2.trojectry[,4]<-((exp(vb2.trojectry[,4]) -1)/(exp(vb2.trojectry[,4]) +1))
# vb1.trojectry[10000,]
# vb2.trojectry[10000,]

vb2.trojectry<-as.data.frame(vb2.trojectry)
vb2.trojectry$Time<-1:15000
vb1.trojectry_melted <- melt(vb2.trojectry, id.vars = "Time")

legend_labels <- c("intercept","log.prop.house",
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
    legend.position = c(0.25, 0.8),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )+ggtitle("HVB-G")+ylim(-0.6,1.9)

pp.con.vbnob




pp=pp.con.vb1+pp.con.vbnob+plot_layout(ncol = 2,nrow = 1)




ggsave(
  "con_vb_elec_MAR.png"
  ,plot = pp,width = 1350,height = 600,units = "px")
