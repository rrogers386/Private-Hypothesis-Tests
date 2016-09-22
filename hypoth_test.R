### Test the DP Hypothesis Tests ###

source("DP_GOF_test_functions.R")
source("DP_independence_test_functions.R")

######################################
############# GOF Test ###############
######################################
n<-1000
d<-4
null_p<-rep(1/d,d)
data<-rmultinom(1,n,null_p)
### Privacy ###
epsilon<-0.1
delta<-10^(-6)

#Gaussian#
noise<-rnorm(d,0,2*sqrt(log(2/delta))/epsilon)
#Laplace#
#noise<-rexp(d,epsilon/2)-rexp(d,epsilon/2)
priv_data<-data+noise


## Return the p_value for this data table ##
# Asymptotic Based Test #
DP_asympt_GOF_test(n,priv_data,epsilon,delta,null_p)



######################################
#### Independence Test ###############
######################################

### Generate private contingency table ###

n<-3000 # Sample Size
r<-3 # rows 
c<-4 # columns
pi1<-rep(1/r,r)
pi2<-rep(1/c,c)
data<-matrix(rmultinom(1,n,as.vector(pi1%o%pi2)),r,c)
### Privacy ###
epsilon<-0.1
delta<-10^(-6)

#Gaussian#
noise<-rnorm(r*c,0,2*sqrt(log(2/delta))/epsilon)
#Laplace#
#noise<-matrix(r*c,rexp(epsilon/2)-rexp(r*c,epsilon/2),r,c)
priv_data<-data+noise

## Return the p_value for this data table ##
# Asymptotic Based Test #
DP_asympt_independence_test(n,priv_data,epsilon,delta)


