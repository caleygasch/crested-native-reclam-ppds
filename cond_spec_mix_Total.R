### Bivariate posterior predictive models demonstration code
### Code prepared by S Huzurbazar, C Gasch, University of Wyoming
### September 14, 2015

### Analysis for Total microbial abundance
### Log scale for C1, C2, C3.
### NOTE: GIVEN THE RANGE OF VALUES, LOG TRANSFORMATION IS ALMOST LINEAR BUT IT DOES TAKE CARE 
### OF NEGATIVE VALUES
### C1 and UD use BUGS...so read in data differently, set directory, transpose later

## Load libraries
library(rjags)
library(R2OpenBUGS)
library(coda)
library(emdbook)

## load the data
basic_data<-read.csv(file="CW_Data.csv",header=T)
## list of variables to be used is Total (Total microbial)

## Define variables used 
attach(basic_data)

### Treatments: 1=Native, 2=Crested, 3=Undisturbed
### Age: for native, ignore Age=1 (1.5 yrs), use Ages =2 and 3
### Depth: 1 for 0-5 and 2 for 5-15

## compute correlations across depths
sink("corr.txt")

W2_Total_1<- Total[Treatment==1 & Age==2 & Depth==1]
W2_Total_2<- Total[Treatment==1 & Age==2 & Depth==2]
W2_Total<-rbind(W2_Total_1,W2_Total_2)
print("W2")
print(cor(cbind(W2_Total_1,W2_Total_2),use="na.or.complete"))

W3_Total_1<- Total[Treatment==1 & Age==3 & Depth==1]
W3_Total_2<- Total[Treatment==1 & Age==3 & Depth==2]
W3_Total<-rbind(W3_Total_1,W3_Total_2)
print("W3")
print(cor(cbind(W3_Total_1,W3_Total_2)))

### correlations are transposed since this site has missing obs and will use BUGS ###
C1_Total_1<- Total[Treatment==2 & Age==1 & Depth==1]
C1_Total_2<- Total[Treatment==2 & Age==1 & Depth==2]
C1_Total<-cbind(C1_Total_1,C1_Total_2)
print("C1")
print(cor(C1_Total,use="na.or.complete"))

C2_Total_1<- Total[Treatment==2 & Age==2 & Depth==1]
C2_Total_2<- Total[Treatment==2 & Age==2 & Depth==2]
C2_Total<-rbind(C2_Total_1,C2_Total_2)
print("C2")
print(cor(cbind(C2_Total_1,C2_Total_2)))

C3_Total_1<- Total[Treatment==2 & Age==3 & Depth==1]
C3_Total_2<- Total[Treatment==2 & Age==3 & Depth==2]
C3_Total<-rbind(C3_Total_1,C3_Total_2)
print("C3")
print(cor(cbind(C3_Total_1,C3_Total_2)))

### correlations are transposed since this set will use BUGS ###
UD_Total_1<- Total[Treatment==3 & Age==1 & Depth==1]
UD_Total_2<- Total[Treatment==3 & Age==1 & Depth==2]
UD_Total<-cbind(UD_Total_1,UD_Total_2)
print("UD")
print(cor(UD_Total,use="na.or.complete"))

sink()

################  Parameters to be set #################

n.iter=15000
n.chains=3
n.adapt=5000
n.thin=2
n.iter.bugs=n.adapt + n.iter/n.thin

################ Models for Total #######################

### W2_Total ; rjags model on original scale ###
Y<-W2_Total
W2_Total_model<- jags.model("cond_bvn_model.txt", data=list("Y"=Y), inits=list('mu'=c(1,1)), 
                            n.chains=n.chains, n.adapt=n.adapt)
W2_Total_samples<-coda.samples(W2_Total_model, c('mu', 'tau','Ypred','sigma','rho'),
                               n.iter=n.iter, n.thin=n.thin,n.adapt=n.adapt)

### W3_Total  ; rjags model on original scale ###
Y<-W3_Total
W3_Total_model<- jags.model("cond_bvn_model.txt", data=list("Y"=Y), inits=list('mu'=c(1,1)),  
                            n.chains=n.chains, n.adapt=n.adapt)
W3_Total_samples<-coda.samples(W3_Total_model, c('mu', 'tau','Ypred','sigma','rho'),
                               n.iter=n.iter, n.thin=n.thin,n.adapt=n.adapt)

### C1_Total ; incomplete data so bugs model, log scale #### 
Y<-log(C1_Total)
model.file <- file.path("cond_bugs_bvn_model_ln.txt")
data <- list ("Y")
parameters <- c("mu", "omega","tau", "Ypredln","Ypred","sigma","rho")
inits<-function(){
  list(mu=rnorm(1,0,10),omega=c(1,1))
}
C1_Total.sim <- bugs(data, inits, parameters, n.iter=n.iter.bugs, model.file,
                     n.chains=n.chains,  n.burnin=n.adapt, n.thin=n.thin, DIC=TRUE, 
                     ## SET WD!!!
                     working.directory =  "./",codaPkg=TRUE)
C1_Total_samples<-read.bugs(c("CODAchain1.txt","CODAchain2.txt","CODAchain3.txt"))

### C2_Total  ; rjags model on log scale ##### 
Y<-log(C2_Total)
C2_Total_model<- jags.model("cond_bvn_model_ln.txt", data=list("Y"=Y), inits=list('mu'=c(0,0)),   
                            n.chains=n.chains, n.adapt=n.adapt)
C2_Total_samples<-coda.samples(C2_Total_model, c('mu', 'tau','Ypredln','Ypred','sigma','rho'),
                               n.iter=n.iter,n.thin=n.thin,n.adapt=n.adapt)

####### C3_Total ; rjags model on log scale #####
Y<-log(C3_Total)
C3_Total_model<- jags.model("cond_bvn_model_ln.txt", data=list("Y"=Y), inits=list('mu'=c(0,0)),   
                            n.chains=n.chains, n.adapt=n.adapt)
C3_Total_samples<-coda.samples(C3_Total_model, c('mu', 'tau','Ypredln','Ypred','sigma','rho'),
                               n.iter=n.iter,n.thin=n.thin,n.adapt=n.adapt)

### UD_Total ; incomplete data so bugs model on original scale ###
Y<-UD_Total
model.file <- file.path("cond_bugs_bvn_model.txt")
data <- list ("Y")
parameters <- c("mu", "tau", "Ypred","sigma","rho")
inits<-function(){
  list(mu=rnorm(2,2,10),Ypred=rnorm(2,0,10),omega=c(1,1))
}
UD_Total.sim <- bugs(data, inits, parameters, n.iter=n.iter.bugs, model.file,
                     n.chains=n.chains,  n.burnin=n.adapt, n.thin=n.thin, DIC=TRUE, 
                     ## SET WD!!!
                     working.directory =  "./", codaPkg=TRUE)
UD_Total_samples<-read.bugs(c("CODAchain1.txt","CODAchain2.txt","CODAchain3.txt"))

##### compile summaries
sink("summary_log_Total.txt")
print("Total on log scale")
print("W2 Total")
print(summary(W2_Total_samples))

print("W3 Total")
print(summary(W3_Total_samples))

print("C1 Total")
print(summary(C1_Total_samples))

print("C2 Total")
print(summary(C2_Total_samples))

print("C3 Total")
print(summary(C3_Total_samples))

print("UD Total")
print(summary(UD_Total_samples))

sink()

######## PLOTS of Predictive densities ##########

par(mfrow=c(3,2))

w2_lines<-HPDregionplot(W2_Total_samples,vars=1:2,prob=0.95)
w3_lines<-HPDregionplot(W3_Total_samples,vars=1:2,prob=0.95)
c1_lines<-HPDregionplot(C1_Total_samples,vars=1:2,prob=0.95)
c2_lines<-HPDregionplot(C2_Total_samples,vars=1:2,prob=0.95)
c3_lines<-HPDregionplot(C3_Total_samples,vars=1:2,prob=0.95)
ud_lines<-HPDregionplot(UD_Total_samples,vars=1:2,prob=0.95)

par(mfrow=c(1,1))

pdf("contours.pdf")
par(fig=c(0,1,.7,1))
par(mar=c(0,5,2,2))

plot(ud_lines[[1]]$x,ud_lines[[1]]$y,xaxt="n", xlim=c(0,50),
ylim=c(0,60),type="l",lty=1,
xlab="Predicted Total at 0-5 cm", 
ylab="nmol/g soil: 5-15 cm",main=""
##,main="Total Microbial Abundance: 95% high posterior predictive regions"
)
text(10,40,substitute("95 % HPD contours"))
text(40,10,substitute("horizontal: 0-5cm"))
lines(w2_lines[[1]]$x,w2_lines[[1]]$y,lty=1,col="2")
lines(w3_lines[[1]]$x,w3_lines[[1]]$y,lty=1,col="3")
lines(c1_lines[[1]]$x,c1_lines[[1]]$y,lty=6,col="4")
lines(c2_lines[[1]]$x,c2_lines[[1]]$y,lty=6,col="9")
lines(c3_lines[[1]]$x,c3_lines[[1]]$y,lty=6,col="6")

#legend("topleft",c("UD","W14","W26","C11","C16","C29"),lty=c(1,1,1,6,6,6),lwd=c(2,2,2,2,2,2),
#       col=c(1,2,3,4,9,6))

#dev.off()

######## PLOTS of Predictive densities ##########

#pdf("Total_plots.pdf")
##par(mfrow=c(2,1))
### upper depth (0-5) $$$

Y1_W2_pred<-unlist(W2_Total_samples[,1,drop=F])
Y1_W3_pred<-unlist(W3_Total_samples[,1,drop=F])
Y1_C1_pred<-unlist(C1_Total_samples[,1,drop=F])
Y1_C2_pred<-unlist(C2_Total_samples[,1,drop=F])
Y1_C3_pred<-unlist(C3_Total_samples[,1,drop=F])
Y1_UD_pred<-unlist(UD_Total_samples[,1,drop=F])

par(fig=c(0,1,.4,.7),new=T)
par(mar=c(0,5,0,2))

plot(density(Y1_UD_pred),lty=1,lwd=2,col=1,xaxt="n",
xlim=c(0,50),
ylim=c(0,max(density(Y1_W2_pred)$y,density(Y1_W3_pred)$y,density(Y1_C1_pred)$y,density(Y1_C2_pred)$y,
             density(Y1_C3_pred)$y,density(Y1_UD_pred)$y)), 
##xlab="nmol/g soil", 
ylab="predictive density",main=""
#,main="(b) Total Microbial Abundance, depth= 0-5"
)
text(10,0.2,substitute("depth = 0-5 cm"))
lines(density(Y1_W2_pred,from=0.0),lty=1,lwd=2, col=2)
lines(density(Y1_W3_pred,from=0.0),lty=1,lwd=2, col=3)
lines(density(Y1_C1_pred,from=0.0),lty=6,lwd=2,col=4)
lines(density(Y1_C2_pred,from=0.0),lty=6,lwd=2,col=9)
lines(density(Y1_C3_pred,from=0.0),lty=6,lwd=2,col=6)
legend("center",c("UD","W14","W26","C11","C16","C29"),lty=c(1,1,1,6,6,6),lwd=c(2,2,2,2,2,2),
       col=c(1,2,3,4,9,6))

##### lower depth #########

Y2_W2_pred<-unlist(W2_Total_samples[,2,drop=F])
Y2_W3_pred<-unlist(W3_Total_samples[,2,drop=F])
Y2_C1_pred<-unlist(C1_Total_samples[,2,drop=F])
Y2_C2_pred<-unlist(C2_Total_samples[,2,drop=F])
Y2_C3_pred<-unlist(C3_Total_samples[,2,drop=F])
Y2_UD_pred<-unlist(UD_Total_samples[,2,drop=F])

par(fig=c(0,1,0,0.4),new=T)
par(mar=c(5,5,0,2))
ymax=max(density(Y2_W2_pred)$y,density(Y2_W3_pred)$y,density(Y2_C1_pred)$y,density(Y2_C2_pred)$y,
         density(Y2_C3_pred)$y,density(Y2_UD_pred)$y)
plot(density(Y2_UD_pred), lty=1,lwd=2,col=1, xlim=c(0,50), ylim=c(0,ymax+.1),xlab="nmol/g soil", 
     ylab="predictive density",main="")
###,main="(c) Total Microbial Abundance, depth=5-15 cm")
text(40,0.5,substitute("depth = 5-15 cm"))
lines(density(Y2_W2_pred,from=0.0),lty=1,lwd=2, col=2)
lines(density(Y2_W3_pred,from=0.0),lty=1,lwd=2, col=3)
lines(density(Y2_C1_pred,from=0.0),lty=6,lwd=2,col=4)
lines(density(Y2_C2_pred,from=0.0),lty=6,lwd=2,col=9)
lines(density(Y2_C3_pred,from=0.0),lty=6,lwd=2,col=6)
#legend("center",c("UD","W14","W26","C11","C16","C29"),lty=c(1,1,1,6,6,6),lwd=c(2,2,2,2,2,2),
#       col=c(1,2,3,4,9,6))

dev.off()

#############################################################################
pdf("log_Total_plots.pdf")
par(mfrow=c(2,1))
### upper depth (0-5) $$$

Y1_W2_pred<-unlist(W2_Total_samples[,1,drop=F])
Y1_W3_pred<-unlist(W3_Total_samples[,1,drop=F])
Y1_C1_pred<-unlist(C1_Total_samples[,1,drop=F])
Y1_C2_pred<-unlist(C2_Total_samples[,1,drop=F])
Y1_C3_pred<-unlist(C3_Total_samples[,1,drop=F])
Y1_UD_pred<-unlist(UD_Total_samples[,1,drop=F])

plot(density(Y1_W2_pred),xlim=c(0,60),
      ylim=c(0,max(density(Y1_W2_pred)$y,density(Y1_W3_pred)$y,density(Y1_C1_pred)$y,
                   density(Y1_C2_pred)$y,density(Y1_C3_pred)$y,density(Y1_UD_pred)$y)),
     xlab="Total", ylab="predictive density",
     main="Predictive distributions for back transformed log Total at 0-5")
lines(density(Y1_W3_pred),lty=2,col=2)
lines(density(Y1_C1_pred),lty=3,col=3)
lines(density(Y1_C2_pred),lty=4,col=4)
lines(density(Y1_C3_pred),lty=5,col=5)
lines(density(Y1_UD_pred),lty=6,col=6)
legend("topright",c("W2","W3","C1","C2","C3","UD"),lty=c(1:6),col=c(1:6))

##### lower depth #########

Y2_W2_pred<-unlist(W2_Total_samples[,2,drop=F])
Y2_W3_pred<-unlist(W3_Total_samples[,2,drop=F])
Y2_C1_pred<-unlist(C1_Total_samples[,2,drop=F])
Y2_C2_pred<-unlist(C2_Total_samples[,2,drop=F])
Y2_C3_pred<-unlist(C3_Total_samples[,2,drop=F])
Y2_UD_pred<-unlist(UD_Total_samples[,2,drop=F])
plot(density(Y2_W2_pred),xlim=c(0,60), ylim=c(0,max(density(Y2_W2_pred)$y,density(Y2_W3_pred)$y,
                                                    density(Y2_C1_pred)$y,density(Y2_C2_pred)$y,
                                                    density(Y2_C3_pred)$y,density(Y2_UD_pred)$y)),
     xlab="Total", ylab="predictive density",
     main="Predictive distributions for back transformed log Total at 5-15")
lines(density(Y2_W3_pred),lty=2,col=2)
lines(density(Y2_C1_pred),lty=3,col=3)
lines(density(Y2_C2_pred),lty=4,col=4)
lines(density(Y2_C3_pred),lt=5,col=5)
lines(density(Y2_UD_pred),lty=6,col=6)
legend("topright",c("W2","W3","C1","C2","C3","UD"),lty=c(1:6),col=c(1:6))

dev.off()

################# hpd intervals #############

sink("log Total_hpd.txt")
print("log Total Weight data")
print("0-5 W2")
print(summary(W2_Total_samples)$quantiles[1,] ) 

print("0-5 W3")
print(summary(W3_Total_samples)$quantiles[1,] )

print("0-5 C1")
print(summary(C1_Total_samples)$quantiles[1,] )

print("0-5 C2")
print(summary(C2_Total_samples)$quantiles[1,] )

print("0-5 C3")
print(summary(C3_Total_samples)$quantiles[1,] )

print("0-5 UD")
print(summary(UD_Total_samples)$quantiles[1,] )

###

print("5-15 W2")
print(summary(W2_Total_samples)$quantiles[2,] )

print("5-15W3")
print(summary(W3_Total_samples)$quantiles[2,] )

print("5-15 C1")
print(summary(C1_Total_samples)$quantiles[2,] )

print("5-15 C2")
print(summary(C2_Total_samples)$quantiles[2,] )

print("5-15 C3")
print(summary(C3_Total_samples)$quantiles[2,] )

print("5-15 UD")
print(summary(UD_Total_samples)$quantiles[2,])

sink()

## end script