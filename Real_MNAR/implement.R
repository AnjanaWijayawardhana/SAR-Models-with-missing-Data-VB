source("source.R")
############################# Data set

# Extract W

hlw <- nb2listw(e80_queen,zero.policy = TRUE)
# Convert listw to a dense matrix
W_dense <- listw2mat(hlw)
# Convert dense matrix to sparse CsparseMatrix
W <- Matrix(W_dense, sparse = TRUE)



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
head(X)

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
head(Xst)
X=Xst

# Make missing values

psi=c(1.4,rep(0.5,1),-0.1) #set psi's (changing these values, we can control the missing value %)

w=W
x=X
y=log.y
y=as.matrix(y)

splitting=my_splitter_all_xm(x=x,xm=as.matrix(x[,1]),y=y,w=w,psi=psi)
m=splitting$m
x=splitting$X
x=as.matrix(x)
yo=splitting$yo
yu=splitting$yu
w=splitting$W
xm=splitting$xm
# x==xm
mean(m) # This is the missing value proportion
nu=sum(m)
nu

# Set starting values, yu.start and theta.start

#-starting values for thetas

n<-ncol(w)
no<-n-nu
dim(x)
length(yo)
reg.model<-estimate_regression(x=x,yo=yo,no=no)
reg.model.summary<-summary(reg.model)
sigma2y0<-(reg.model.summary$sigma)^2
beta0<-reg.model.summary$coefficients

prob1<-sum(m)/length(m)

psi0_0<-log(prob1/(1-prob1))
# paras.st<-c(as.numeric(beta0[,1]),sigma2y0,0.001,psi0_0,rep(0.1,2))
paras.st<-c(as.numeric(beta0[,1]),sigma2y0,0.01,psi0_0,rep(0.1,ncol(xm)+1))
wpluswt=(w+t(w))
wtw=t(w)%*%(w)

#-starting values for yu

N1=50
bsize=100 # set block size
k=floor(nu/bsize) # # of full blocks
reminder<-nu%%bsize # length of remaining incomplete block
reminder
blocks=k+ifelse(reminder>0,1,0)
listw_x_m<-make.w_x_(w=w,x=cbind(1,x),xm=cbind(1,xm),m=m,k=k,reminder =reminder)
startyu<-rep(0.1,nu)

tic()
mysamples<-sim_yugyotheta_(yo=yo,listw_x_m=listw_x_m,
                           yu.start=startyu,rc = ncol(x)+1,
                           paras=paras.st,N1=N1,k=k,bsize=bsize,reminder=reminder)

toc()

(mysamples$accepts)/(N1*blocks)
dim((mysamples$yugyom.chain))
yu.start<-(mysamples$yugyom.chain)[N1,] # We give the starting value as yu,
plot(yu,yu.start)


########################### Implement VB algorithms ########################################


N=15000

## HVB-AllB
bsize=floor(nu*0.1) # selecting block size to be nu*0.25
bsize.upall<-bsize

tic()
fit_thetaAug_update_all<-MNAR_VB_thetaAug_fullB_updating(x=x,xm=xm,m=m,
                                                         w=w,yo=yo,startyu=yu.start,p=4,
                                                         start_theta=paras.st,N=N,N1=20,bsize=bsize.upall)
toc()


### HVB-3B
# N=10000

bsize.updateR3<-bsize
tic()
fit_thetaAug_updateR3<-MNAR_VB_thetaAug_R3B_updating(x=x,xm=xm,m=m,
                                                     w=w,yo=yo,startyu=yu.start,p=4,
                                                     start_theta=paras.st,N=N,N1=20,
                                                     bsize=bsize.updateR3) 
toc()



## JVB

st=c(yu.start,paras.st)


tic()
fit_thetayu<-MNAR_thetayu(x=x,xm=xm,m=m,w=w,yo=yo,p=4,
                          start_theta=st,N=N)
toc()

rnames=c("intercept",colnames(X),"sigma2y","rho","psi0","psix","psiy")

###### Sample from posteriors of theta

# JVB

sample.size<-10000
nparas<-length(fit_thetaAug_update_all$mu)
yu.theta.sample<-matrix(rep(0,sample.size*(nparas+nu)),ncol = nparas+nu)


for(i in 1:sample.size){
  yu.theta.sample[i,]<-simulate_VB_yu.theta(mu=fit_thetayu$mu,x=x,
                                            B=fit_thetayu$B,
                                            d=fit_thetayu$d)
  
}

theta.sample.JVB<-yu.theta.sample[,-(1:nu)]

colnames(theta.sample.JVB)=rnames
#post mean and sd
round(apply(theta.sample.JVB, 2, mean),4)
round(apply(theta.sample.JVB, 2, sd),4)


# HVB-AllB

theta.sample.upall<-matrix(rep(0,sample.size*(nparas)),ncol = nparas)

for(i in 1:sample.size){
  theta.sample.upall[i,]<-simulate_VB_theta(mu=fit_thetaAug_update_all$mu,
                                            B=fit_thetaAug_update_all$B,
                                            d=fit_thetaAug_update_all$d,
                                            x=x)
  
}

colnames(theta.sample.upall)=rnames
#post mean and sd
round(apply(theta.sample.upall, 2, mean),4)
round(apply(theta.sample.upall, 2, sd),4)



#HVB-3B

theta.sample.up3<-matrix(rep(0,sample.size*(nparas)),ncol = nparas)

for(i in 1:sample.size){
  theta.sample.up3[i,]<-simulate_VB_theta(mu=fit_thetaAug_updateR3$mu,
                                          B=fit_thetaAug_updateR3$B,
                                          d=fit_thetaAug_updateR3$d,
                                          x=x)
  
}


colnames(theta.sample.up3)=rnames
#post mean and sd
round(apply(theta.sample.up3, 2, mean),4)
round(apply(theta.sample.up3, 2, sd),4)


############## simulate and plot Missing values


#Simulate HVB missing values


N1=20
bsize=20
k=floor(nu/bsize) # # of full blocks
reminder<-nu%%bsize # length of remaining incomplte block
blocks=k+ifelse(reminder>0,1,0)
listw_x_m<-make.w_x_(w=w,x=cbind(1,x),xm=cbind(1,xm),m=m,k=k,reminder =reminder)


tic()
yu_Aug_all<-sim_yugyotheta_(yo=yo,listw_x_m=listw_x_m,
                            yu.start=startyu,rc = ncol(x)+1,
                            paras=apply(theta.sample.upall,2,mean),N1=N1,k=k,bsize=bsize,reminder=reminder)

toc()




tic()
yu_Aug_3B<-sim_yugyotheta_(yo=yo,listw_x_m=listw_x_m,
                           yu.start=startyu,rc = ncol(x)+1,
                           paras=apply(theta.sample.up3,2,mean),N1=N1,k=k,bsize=bsize,reminder=reminder)

toc()




# VB means vs true yu 

legend_labels<-c("JVB","HVB-AllB","HVB-3B"
)

data <- data.frame(
  x=yu,  # Replace with your x values
  vb1=apply(yu.theta.sample[,1:nu],2,mean),  # Replace with your y values
  vb2allup=apply((yu_Aug_all$yugyom.chain)[(N1/2):N1,],2,mean),
  vb23B=apply((yu_Aug_3B$yugyom.chain)[(N1/2):N1,],2,mean)
)


pp1 <- ggplot(data, aes(x = x)) +
  geom_point(aes(y = vb1, color = "vb1"), size = 2) +
  geom_point(aes(y = vb2allup, color = "vb2allup"), size = 2) +
  geom_point(aes(y = vb23B, color = "vb23B"), size = 2) + 
  labs(
    title = "",
    x = "True missing value",
    y = "VB posterior means"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue","red",
                                "purple"), labels = legend_labels)  +
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    legend.position = c(0.25,0.8),  # Position the legend at the top of the plot
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a border around the plot
    
    axis.text = element_text(size = 22),  # Adjust the size as needed
    axis.title = element_text(size = 22),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5, size = 25),  # Adjust the size as needed
    legend.text = element_text(size = 22)  # Adjust the size as needed
  ) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "green") +
  guides(color = guide_legend(title = ""))  # legend title

pp1  # Display the plot

# calculate MSEs
sum((yu-apply((yu_Aug_all$yugyom.chain)[(N1/2):N1,],2,mean))^2) #HVB-AllB
sum((yu-apply((yu_Aug_3B$yugyom.chain)[(N1/2):N1,],2,mean))^2) #HVB-3B
sum((yu-apply(yu.theta.sample[,1:nu],2,mean))^2) #JVB


#############  Generate density Plots



name.list<-c("VB-yutheta","VB-upall","VB-up3")

values=c("VB-yutheta"="blue","VB-upall"="red","VB-up3"="purple")





df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(yu.theta.sample[,(nu+1)],theta.sample.upall[,1],
            theta.sample.up3[,1])
)


pp1 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5) +
  theme_minimal() +
  scale_x_continuous(limits = c(-0.75, -0.5)) +
  scale_y_continuous(limits = c(0, 51)) +
  scale_color_manual(values = values, breaks = c("VB-yutheta","VB-upall","VB-up3"),
                     labels = c("JVB","HVB-AllB","HVB-3B")) +  # Specify breaks
  labs(
    x = "intercept",
    y = "Density",
    color = "") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.25, 0.82),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 22),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    legend.text.align = 0# Adjust the size as needed
  )  # Adjust the size as needed

pp1

head(X)
#log.house, n_baths

df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(yu.theta.sample[,(nu+3)],theta.sample.upall[,3],theta.sample.up3[,3])
)

ppblog.house <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(0.1, 0.7)) +
  scale_y_continuous(limits = c(0, 12)) +
  scale_color_manual(values = values, breaks = c("VB-yutheta","VB-upall","VB-up3"),
                     labels = c("JVB","HVB-AllB","HVB-3B"))+
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
  )  # Adjust the size as needed

ppblog.house




#sigma2y
nbeta<-ncol(X)+1

df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(yu.theta.sample[,(nu+nbeta+1)],theta.sample.upall[,nbeta+1],
            theta.sample.up3[,nbeta+1])
)

pp2 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 0.035)) +
  scale_y_continuous(limits = c(0, 750)) +
  scale_color_manual(values = values, breaks = c("VB-yutheta","VB-upall","VB-up3"),
                     labels = c("JVB","HVB-AllB","HVB-3B")) +
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
    legend.text = element_text(size = 5))  # Adjust the size as needed

pp2


#rho

df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(yu.theta.sample[,(nu+nbeta+2)],theta.sample.upall[,nbeta+2],
            theta.sample.up3[,nbeta+2])
)

pp3 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 25)) +
  scale_color_manual(values = values, breaks = c("VB-yutheta","VB-upall","VB-up3"),
                     labels = c("JVB","HVB-AllB","HVB-3B"))+
  labs(
    x = expression(rho),
    y = "Density",
    color = "") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.5, 0.72),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 22),  # Adjust the size as needed
    legend.text = element_text(size = 5))  # Adjust the size as needed

pp3


#psi0
df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c( yu.theta.sample[,(nu+nbeta+3)],
             theta.sample.upall[,nbeta+3],theta.sample.up3[,nbeta+3])
)

pp4 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 3)) +
  scale_y_continuous(limits = c(0, 2.8)) +
  scale_color_manual(values = values, breaks = c("VB-yutheta","VB-upall","VB-up3"),
                     labels = c("JVB","HVB-AllB","HVB-3B")) +
  labs(
    x = expression(psi[0]),
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
    legend.text = element_text(size = 5))  # Adjust the size as needed

pp4

#psiy


df <- data.frame(
  Group = rep(name.list, each = sample.size),
  Value = c(yu.theta.sample[,(nu+nbeta+5)],
            theta.sample.upall[,nbeta+5],theta.sample.up3[,nbeta+5])
)


pp5 <- ggplot(df, aes(x = Value, color = Group)) +
  geom_density(fill = NA, size = 0.5,show.legend = FALSE) +
  theme_minimal() +
  scale_x_continuous(limits = c(-1, 1.6)) +
  scale_y_continuous(limits = c(0, 1.8)) +
  scale_color_manual(values = values, breaks = c("VB-yutheta","VB-upall","VB-up3"),
                     labels = c("JVB","HVB-AllB","HVB-3B")) +
  labs(
    x = expression(psi[y]),
    y = "Density",
    color = "") +
  theme(
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.12, 0.72),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 22),  # Adjust the size as needed
    legend.text = element_text(size = 5))  # Adjust the size as needed

pp5


pp=pp1+ppblog.house+pp2+pp3+pp4+pp5+plot_layout(ncol = 3,nrow = 2)
ggsave(file = "kernal_SEM_MNAR_elec.png",
       pp,width = 1350,height = 768,units = "px") # If we set these, dimensions, the hwole
# frames as well as individul graps make similar sizes



##### Convergence plots 

## JVB lower bound


# Convert time series to dataframe
ts_data <- data.frame(Iteration = time(fit_thetayu$L.bound), 
                      Value = fit_thetayu$L.bound)

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

# HVB-AllB

# all
vb2.trojectry=(fit_thetaAug_update_all$all_paras)[,c(1,9,10,11)]
vb2.trojectry[,2]<-exp(vb2.trojectry[,2])
vb2.trojectry[,3]<-((exp(vb2.trojectry[,3]) -1)/(exp(vb2.trojectry[,3]) +1))
# vb1.trojectry[10000,]


vb2.trojectry<-as.data.frame(vb2.trojectry)
vb2.trojectry$Time<-1:15000
vb1.trojectry_melted <- melt(vb2.trojectry, id.vars = "Time")

legend_labels <- c("intercept",
                   expression(sigma[e]^2), expression(rho), expression(psi[0]))

# Define line types
line_types <- c("solid", "dashed", "dotted", "dotdash")

pp.con.vball <- ggplot(data = vb1.trojectry_melted, aes(x = Time, y = value, color = variable, linetype = variable)) +
  geom_line(size = 0.25) +
  scale_color_manual(values = rep("red", length(legend_labels)), name = NULL, labels = legend_labels) +
  scale_linetype_manual(values = line_types, name = NULL, labels = legend_labels) + # Specify line types and remove linetype legend title
  labs(title = "",
       x = "Iteration", y = "Value") +
  theme_minimal() +
  theme(
    
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.4, 0.8),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    # legend.direction = "horizontal",
    legend.text.align = 0# Adjust the size as needed
  )+ggtitle("HVB-AllB")+ylim(-1,4)

pp.con.vball



# 3B


vb2.trojectry=(fit_thetaAug_updateR3$all_paras)[,c(1,9,10,11)]
vb2.trojectry[,2]<-exp(vb2.trojectry[,2])
vb2.trojectry[,3]<-((exp(vb2.trojectry[,3]) -1)/(exp(vb2.trojectry[,3]) +1))


vb2.trojectry<-as.data.frame(vb2.trojectry)
vb2.trojectry$Time<-1:15000
vb1.trojectry_melted <- melt(vb2.trojectry, id.vars = "Time")



# Define line types
line_types <- c("solid", "dashed", "dotted", "dotdash")

pp.con.vb3 <- ggplot(data = vb1.trojectry_melted, aes(x = Time, y = value, color = variable, linetype = variable)) +
  geom_line(size = 0.25) +
  scale_color_manual(values = rep("purple", length(legend_labels)), name = NULL, labels = legend_labels) +
  scale_linetype_manual(values = line_types, name = NULL, labels = legend_labels) + # Specify line types and remove linetype legend title
  labs(title = "",
       x = "Iteration", y = "Value") +
  theme_minimal() +
  theme(
    
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.4, 0.8),
    panel.border = element_blank(),
    panel.background = element_rect(color = "black", fill = NA, size = 0.5),
    
    axis.text = element_text(size = 5),  # Adjust the size as needed
    axis.title = element_text(size = 5),  # Adjust the size as needed
    plot.title = element_text(hjust = 0.5,size = 5),  # Adjust the size as needed
    legend.text = element_text(size = 5),  # Adjust the size as needed
    legend.key.size = unit(0.4, "lines"),
    # legend.direction = "horizontal",
    legend.text.align = 0# Adjust the size as needed
  )+ggtitle("HVB-3B")+ylim(-1,4)

pp.con.vb3

pp=pp.con.vb1+pp.con.vball+pp.con.vb3+plot_layout(ncol = 2,nrow = 2)


ggsave(
  "con_vb_MNAR_elec.png"
  ,plot = pp,width = 1350,height = 1100,units = "px")
nu




