# --------------------------------------------------------------------------------------------------
#	
#	Robert-Casella - Section 7.2, normal example 7.3
#
# --------------------------------------------------------------------------------------------------
require(mcsm)
data(Energy)
# --------------------------------------------------------------------------------------------------
x<-Energy$Girls
n<-length(x)
xbar<-mean(x)
# --------------------------------------------------------------------------------------------------
#
#	Parameters for prior pdf's according to Roberts-Casella
#
# --------------------------------------------------------------------------------------------------
a<-3
b<-3
theta0<-5
tau2<-10
# --------------------------------------------------------------------------------------------------
#
#	Empirical quantities for the given sample
#
# --------------------------------------------------------------------------------------------------
s2tilde<-var(x)
stilde<-sd(x)
round(xbar,3)
#[1] 870.5
round(s2tilde,3)
#[1] 143093.1
round(s2tilde/n,3)
#[1] 8943.317
lambda1<-(s2tilde/n)/(s2tilde/n+tau2)
round(lambda1,3)
#[1] 0.999
# --------------------------------------------------------------------------------------------------
#
#	Empirical quantities for the full population
#
# --------------------------------------------------------------------------------------------------
Energy.Girls<-read.table("Energy.girls.txt",header=TRUE)
X<-Energy.Girls$Girls
N<-length(X)
Xbar<-mean(X)
Sqtilde<-var(X)
round(Xbar,3)
#[1] 917.855
round(Sqtilde,3)
#[1] 143876.3
round(Sqtilde/N,3)
#[1] 1307.966
# --------------------------------------------------------------------------------------------------
#
#	Prior parameters for theta
#
# --------------------------------------------------------------------------------------------------
tau2<-Sqtilde/N
theta0<-Xbar
# --------------------------------------------------------------------------------------------------
#
#	Prior quantities for sigma2: We fit b from the population unit variance
#
# --------------------------------------------------------------------------------------------------
a<-3
#Prior.sigma2.expectation<-b/(a-1)
Prior.sigma2.expectation<-Sqtilde/N
b<-Prior.sigma2.expectation*(a-1)
round(b,3)
#[1] 2615.932
Prior.sigma2.variance<-b^{2}/((a-1)^{2}*(a-2))
Prior.sigma2.variance<-b^{2}/((a-1)^{2}*(a-2))
round(Prior.sigma2.variance,3)
#[1] 1710775
# --------------------------------------------------------------------------------------------------
#
#	Number of samples to generate
#
# --------------------------------------------------------------------------------------------------
nsim<-10^4
# --------------------------------------------------------------------------------------------------
#
#	Run the Gibbs sampler
#
# --------------------------------------------------------------------------------------------------
sh1<-(n/2)+a
theta<-rep(0,nsim)		                #init arrays
sigma2<-rep(0,nsim)      		        #init arrays
sigma2[1]<-1/rgamma(1,shape=a,rate=b)   #init chains
B<-sigma2[1]/(sigma2[1]+n*tau2)
theta[1]<-rnorm(1,m=B*theta0+(1-B)*xbar,sd=sqrt(tau2*B))
for (i in 2:nsim){
   B<-sigma2[i-1]/(sigma2[i-1]+n*tau2)
   theta[i]<-rnorm(1,m=B*theta0+(1-B)*xbar,sd=sqrt(tau2*B))
   ra1<-(1/2)*(sum((x-theta[i])^2))+b
   sigma2[i]<-1/rgamma(1,shape=sh1,rate=ra1)
   }
# --------------------------------------------------------------------------------------------------
#
#	Summarize and plot results
#
# --------------------------------------------------------------------------------------------------
summary(theta)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  784.1   887.8   909.7   909.9   932.1  1029.0 
summary(sigma2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1090   84000  103000  109600  126600  585200 
# --------------------------------------------------------------------------------------------------
log.theta<-log(theta)
log.sigma2<-log(sigma2)
summary(log.theta)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  6.665   6.789   6.813   6.813   6.837   6.936 
summary(log.sigma2)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  6.993  11.340  11.540  11.560  11.750  13.280 
# --------------------------------------------------------------------------------------------------
old.par<-par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(log.theta,nclass=140,col="grey",freq=FALSE,main="",ylab="",xlab=expression(theta))
hist(log.sigma2,nclass=150,col="sienna",freq=FALSE,xlab=expression(sigma^2),main="",ylab="")
par(old.par)
# --------------------------------------------------------------------------------------------------
old.par<-par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(theta,nclass=140,col="grey",freq=FALSE,main="",ylab="",xlab=expression(theta))
hist(sigma2,nclass=150,col="sienna",freq=FALSE,xlab=expression(sigma^2),main="",ylab="")
par(old.par)
# --------------------------------------------------------------------------------------------------

