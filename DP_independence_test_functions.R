#### DP Independence Tests ###
#install.packages("limSolve")
#install.packages("numDeriv")
#install.packages("Matrix")
#install.packages("CompQuadForm")

library("limSolve")
library("numDeriv")
library("Matrix")
library("CompQuadForm")

#### MLE CALCULATIONS #####

# Gamma parameter should only be set >0 if using Laplace noise #
privMLE<-function(n,priv_data,gamma=0){#Data should be a table
  k<-nrow(priv_data)
  ell<-ncol(priv_data)
  priv_vect<-as.vector(priv_data)
  #Calculate MLE for the probability vector.  We do this in two steps:
  if(gamma==0){
    # Step 1, find most likely contingency table that sums to N given the counts have normal noise added to them
    MLE_CT<-matrix(lsei(A=diag(k*ell),B=priv_vect,E=t(rep(1,k*ell)),F=n, G = diag(k*ell),H = rep(0,k*ell))$X,k,ell)
  }else{
    # Step 1, find most likely contingency table that sums to N given the counts have normal noise added to them
    opt_funct<-function(x){
      loss1<-sum(abs(priv_vect-x))
      loss2<-sum((priv_vect-x)^2)
      return((1-gamma)*loss1+gamma*loss2)
    }
    opt_val<-solnp(rep(n/(k*ell),k*ell), #starting values (obviously need to be positive)
                   opt_funct, #function to optimise
                   eqfun=sum, #equality function 
                   eqB=n,   #the equality constraint
                   LB=rep(0,k*ell), #lower bound for parameters i.e. greater than zero
                   UB=rep(n,k*ell))$pars
    opt_table<-matrix(opt_val,k,ell)
  }
  # Step 2, given most likely contingency table, what is the MLE for the probability vector that generated this table
  # given that the variables are independent.
  if(sum(MLE_CT<5)>1){ # In case the MLE table gives some entry with too small of a count #
    MLE_pis<-0
  }else{
    expected_ct<- chisq.test(MLE_CT)$expected# pi1,pi2
    MLE_pis<-list(apply(expected_ct,1,sum)/n,apply(expected_ct,2,sum)/n) #pi1 is [[1]], pi2 is [[2]]
  }
  return(MLE_pis)
}


DP_asympt_independence_test<-function(n,priv_data,epsilon,delta,gamma=0){
  # Chi squared statistic
  Q<-chisq.test(priv_data)$statistic
  
  r<-nrow(priv_data)
  c<-ncol(priv_data)
  ### We add gaussian noise to the data with the following stnd dev. ###
  
  sigma<-2*sqrt(log(2/delta))/epsilon
  
  f<-function(z){
    x<-c(z[1:r-1],1-sum(z[1:r-1]))
    y<-c(z[(r):(r+c-2)],1-sum(z[(r):(r+c-2)]))
    total_matrix<-x%o%y
    total_vect<-as.vector(total_matrix)
    return(total_vect)
  }
  MLE_pi1pi2<-privMLE(n,priv_data,gamma)
  
  pi1<-MLE_pi1pi2[[1]]
  pi2<-MLE_pi1pi2[[2]]
  
  p_vect<-as.vector(pi1%o%pi2)
  
  ind_matrix<-diag(1/sqrt(p_vect))%*%jacobian(f,c(pi1[1:(r-1)],pi2[1:(c-1)]))
  additional_matrix<-ind_matrix%*%solve(t(ind_matrix)%*%ind_matrix)%*%t(ind_matrix)
  covariance_ind<-diag(r*c) - sqrt(p_vect)%*%t(sqrt(p_vect)) - additional_matrix
  
  ### The variable Big_Sigma is the covariance matrix for the vector of random variables where 
  ### the first len components are (X-Np_vect)/sqrt(Np_vect(1-p_vect)) and the second len components 
  ### are the i.i.d. normal random variables ~ N(0,sigma^2) we add for privacy.  We call the 2*len 
  ### dimensional random vector Z.  
  Big_Sigma_ind<-diag(2*(r*c))
  Big_Sigma_ind[1:(r*c),1:(r*c)]<-covariance_ind
  ### We now need to construct the matrix A such that the private chi-squared statistic Q^2_epsilon = W^T A W.
  LT<-diag(r*c)
  RT<-diag(sigma/sqrt(n*p_vect))
  RB<-diag(sigma^2/(n*p_vect))
  A<-rbind(cbind(LT,RT),cbind(RT,RB))
  
  ### Find threshold using quadratic form of normals ###
  # First compute vector of eigenvectors for our singular covariance matrix #
  # Must find B, such that B B^T = Big_Sigma and B is of rank 2*len-1
  rank_val<-rankMatrix(Big_Sigma_ind)[1]
  if(eigen(Big_Sigma_ind)$values[rank_val]<=0){
    rank_val<-rank_val-1
  }
  B_ind<-sapply(1:rank_val,function(i) sqrt(eigen(Big_Sigma_ind)$values[i]) *eigen(Big_Sigma_ind)$vectors[,i])
  lambda_ind<-eigen(t(B_ind)%*%A%*%B_ind)$values
  rank_min<-rankMatrix(A)[1]
  lambda_ind<-lambda_ind[1:rank_min]
  # Numerically find the threshold
  p_val<-imhof(Q,lambda=lambda_ind)$Qq
  return(p_val)
}

