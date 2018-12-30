# --------------------------------------------------------------------------------------------------
#	
#	Robert-Casella - Section 7.2, Beta-binomial
#
# --------------------------------------------------------------------------------------------------
#
#	Beta-binomial joint pdf
#
# --------------------------------------------------------------------------------------------------
betabi=function(x,a,b,n){
   beta(x+a,n-x+b)/(beta(a,b)*beta(x+1,n-x+1)*(n+1))}
# --------------------------------------------------------------------------------------------------
#
#	Gibbs sampler for the Beta-binomial model
#
# --------------------------------------------------------------------------------------------------
nsim<-10^4
n<-15
a<-3
b<-7
X<-rep(0,nsim)
T<-rep(0,nsim)
T[1]<-rbeta(1,a,b)               #initialize the chain
X[1]<-rbinom(1,n,T[1])           #initialize the chain
for (i in 2:nsim){
    X[i]<-rbinom(1,n,T[i-1])
    T[i]<-rbeta(1,a+X[i],n-X[i]+b)
	}
par(mfrow=c(1,2),mar=c(4,4,2,1))
hist(X[2000:nsim],nclass=16,col="grey",freq=FALSE, xlim=c(0,15),main="",xlab="X")
curve(betabi(x,a,b,n),from=0, to=15,col="gold4",lwd=2,add=TRUE)
hist(T[2000:nsim],nclass=134,col="grey",freq=FALSE,xlim=c(0,.8),main="", xlab=expression(theta))
curve(dbeta(x,shape1=a,shape2=b),from=0, to=.8,col="sienna",lwd=2,add=TRUE)