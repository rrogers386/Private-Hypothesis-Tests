#### DP GOF Tests ###
#install.packages("limSolve")
#install.packages("numDeriv")
#install.packages("Matrix")
#install.packages("CompQuadForm")

library("limSolve")
library("numDeriv")
library("Matrix")
library("CompQuadForm")


mychisq<-function(val,n,p){
  return(sum((val - n*p)^2/(n*p)))
}


DP_asympt_GOF_test<-function(n,priv_data,epsilon,delta,null_p){
  Q<-mychisq(priv_data,n,null_p)
  
  ######################################################################################################################
  ### As evidenced in past results, it seems like testing proportions with a chi-squared statistic varies a lot when   ###
  ### noise is added to the vector of counts.  In this script we modify the threshold value for when we have a 1-alpha ###
  ### level of significance.  This involves a quadratic form of mulivariate normals.
  
  ### Hypothesis test H_0: p= p_vect
  len<-length(null_p)
  df<-len-1
  ### We add gaussian noise to the data
  sigma<-2*sqrt(log(2/delta))/epsilon
  ### Our data will be sampled from a mulinomial distribution with parameters N and p_vect
  ### The variable covariance gives the covariance matrix of the vector 
  ### of random variables (X-Np_vect)/sqrt(Np_vect(1-p_vect)) where X ~ Multinomial(N,p_vect)
  covariance<-diag(len) - sqrt(null_p)%*%t(sqrt(null_p))
  
  ### The variable Big_Sigma is the covariance matrix for the vector of random variables where 
  ### the first len components are (X-Np_vect)/sqrt(Np_vect(1-p_vect)) and the second len components 
  ### are the i.i.d. normal random variables ~ N(0,sigma^2) we add for privacy.  We call the 2*len 
  ### dimensional random vector Z.  
  Big_Sigma<-diag(2*len)
  Big_Sigma[1:len,1:len]<-covariance
  ### We now need to construct the matrix A such that the private chi-squared statistic Q^2_epsilon = W^T A W.
  LT<-diag(len)
  RT<-diag(sigma/sqrt(n*null_p))
  RB<-diag(sigma^2/(n*null_p))
  A<-rbind(cbind(LT,RT),cbind(RT,RB))
  
  ### Find threshold using quadratic form of normals ###
  # First compute vector of eigenvectors for our singular covariance matrix #
  # Must find B, such that B B^T = Big_Sigma and B is of rank 2*len-1
  r<-rankMatrix(Big_Sigma)[1]
  B<-sapply(1:r,function(x) sqrt(eigen(Big_Sigma)$values[x]) *eigen(Big_Sigma)$vectors[,x])
  lambda<-eigen(t(B)%*%A%*%B)$values
  r_min<-rankMatrix(A)[1]
  lambda<-lambda[1:r_min]
  # Numerically find the threshold
  p_val<-imhof(Q,lambda=lambda)$Qq
  return(p_val)
}


