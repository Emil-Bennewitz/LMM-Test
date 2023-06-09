

### cdf of the weighted sum of chisquare random variables

#Note: total distribution is simulated "simulations" times <-increasing load with nr factor levels/ranefs 
#gets too much when hundreds of factor levels

#if specify return_data=TRUE, (or something else), then it returns only the simulated data
sum.chi.squares.distr<-function(list.weights,quantile,simulations=10**6,return_data=FALSE){
  if(return_data==FALSE){
    sum=rep(0,simulations)
    for (weight in list.weights){
      sum=sum+weight*rchisq(simulations,1)
    }
    result<-ecdf(sum)(quantile)
    return(result)
  } else{ return(sum) }
}
sum.chi.squares.distr(c(48),100)
sum.chi.squares.distr(rep(1,48),100)

#t<-sum.chi.squares.distr(list(1,1),1)
#print(t)
#pchisq(1,1)

###Auxiliary function

#takes the submatrix consisting of columns pos.start to pos.end and moves it to the end of the matrix, returning modif matrix. 
swap.column<-function(matrix,pos.start,pos.end){
  dim<-dim(matrix)
  n<-dim[1]
  p<-dim[2]
  
  if (pos.start<1 | pos.start>p){
    print("Invalid starting position in swap.block!")
  }
  if (pos.end<pos.start | pos.end>p){
    print("Invalid end position in swap.block!")
  }
  
  if (pos.start==1 && pos.end==p){
    Z.new=matrix
  } else if (pos.start==1 && pos.end<p){
    Z2<-matrix[,pos.start:pos.end]
    Z3<-matrix[,(pos.end+1):p]
    Z.new<-cbind(Z3,Z2)
  } else if (pos.start>1 && pos.end==p){
    #Z1<-matrix[,1:(pos.start-1)]
    #Z2<-matrix[,pos.start:pos.end]
    #Z.new<-cbind(Z2,Z1)
    Z.new<-matrix # I made a mistake here
  } else if (pos.start>1 && pos.end<p){
    Z1<-matrix[,1:(pos.start-1)]
    Z2<-matrix[,pos.start:pos.end]
    Z3<-matrix[,(pos.end+1):p]
    Z.new<-cbind(Z1,Z3,Z2)
  } else {
    print("Something went wrong in swap.block()!")
    browser()
  }
  return(Z.new)
}

#indicated block gets cut pasted to be the last block
#returns sparse matrix, accepts sparse and normal matrix input
swap.block<-function(block.matrix, pos.start, pos.end){
  dim<-dim(block.matrix)
  n=dim[1]
  m=dim[2]
  if (n!=m){
    print("Error in swap.block : matrix must be square!")
  }
  if (pos.start<1 | pos.start>n){
    print("Invalid starting position in swap.block!")
  }
  if (pos.end<pos.start | pos.end>n){
    print("Invalid end position in swap.block!")
  }
  
  if (pos.start==1 && pos.end==n){
    M.new=block.matrix
  } else if (pos.start==1 && pos.end<n){
    M2<-block.matrix[pos.start:pos.end,pos.start:pos.end]
    M3<-block.matrix[(pos.end+1):n,(pos.end+1):n]
    M.new<-bdiag(M3,M2)
  } else if (pos.start>1 && pos.end==n){
    M1<-block.matrix[1:(pos.start-1),1:(pos.start-1)]
    M2<-block.matrix[pos.start:pos.end,pos.start:pos.end]
    M.new<-bdiag(M2,M1)
  } else if (pos.start>1 && pos.end<n){
    M1<-block.matrix[1:(pos.start-1),1:(pos.start-1)]
    M2<-block.matrix[pos.start:pos.end,pos.start:pos.end]
    M3<-block.matrix[(pos.end+1):n,(pos.end+1):n]
    M.new<-bdiag(M1,M3,M2)
  } else {
    print("Something went wrong in swap.block()!")
    browser()
  }
  
  return(M.new)
}


# lmmtest--------------------------------------------------

#square root for symmetric matrix A
square_root<-function(A){
  e<-eigen(A)
  eigenvalues<-e$values
  U<-e$vectors
  n<-dim(A)[[1]]
  
  sqrt_eigenvalues<-sqrt(eigenvalues)
  B<-matrix(rep(0,n^2),nrow=n)
  counter=1
  for (value in sqrt_eigenvalues){
    B<-B+value*U[,counter]%*%t(U[,counter])
    counter<-counter+1
  }
  return(B)
}

pseudo_inverse<-function(A,tol=0){
  e<-eigen(A)
  eigenvalues<-e$values
  U<-e$vectors
  n<-dim(A)[[1]]
  
  B<-matrix(rep(0,n^2),nrow=n)
  counter=1
  for (value in eigenvalues){
    if (value>tol){
      B<-B+(1/value)*U[,counter]%*%t(U[,counter])
    }
    
    counter<-counter+1
  }
  return(B)  
}

#A<-matrix(data=c(1,2,0,2,5,0,0,0,0),byrow=TRUE,nrow=3,ncol=3)
#square_root(A)%*%square_root(A)
#A%*%pseudo_inverse(A)




#fixME: weights are currently ignored for simplicity

#' A simple test for random effects in regression models, LMM case.
#' @description An implementation of the test proposed by Simon Wood in the Biometrika paper "A simple test for random effects in regression models". Link: https://doi.org/10.1093/biomet/ast038
#'
#'
#' @param lmmfit A fitted object from the lmer() function.
#' @param Ztlist.start The first index (inclusive) in Ztlist corresponding to the random effect to be tested.
#' @param Ztlist.end The last index (inclusive) in Ztlist corresponding to the random effect to be tested. For example, giving Ztlist.start=1 and Ztlist.end=2 tells lmmtest that the random effect to be tested consists of elements 1 and 2 in Ztlist. Ztlist is given by getME(lmmfit, "Ztlist"), see examples.
#' @param simulations The number of simulations used to obtain the distribution function of the test statistic. Default is 10^6.
#'
#' @return The p-value of the test. Numeric, can be exactly 0. As simulations determine the test statistic distribution repeatedly every time, results may vary slightly with each execution.
#' @export
#' @import designr Matrix lme4
#' @rdname lmmtest
#' @examples
#' 
#' #Example: LMM, one factor, 100 levels, 5 observations per factor.
#' #Random intercept, random slope, y=3+alpha+z*beta+epsilon, as follows:
#' library(lme4vartest)
#' 
#' set.seed(123)
#' 
#' #create dataframe
#' number.factors<-100 
#' f<-1:number.factors
#' f<-as.factor(f)
#' z1<-rnorm(5*number.factors)
#' intercept<-3
#' test.df<-data.frame(f=f,
#'                     z1=z1)
#' 
#' #simulate responses: y has a strong random slope effect, y_ a very weak effect
#' y<-simLMM(formula=~(1|f)+(0+z1|f),Fixef = c(intercept),VC_sd=list(c(5),c(5),c(2)),data=test.df)
#' y_<-simLMM(formula=~(1|f)+(0+z1|f),Fixef = c(intercept),VC_sd=list(c(5),c(0.4),c(2)),data=test.df)
#' 
#' test.df$y<-y
#' test.df$y_<-y_
#' 
#' #fit models to the simulated data
#' test.model<-lmer(y~(1|f)+(0+z1|f),data=test.df)
#' 
#' 
#' #looking at this object, we can see the different components of the Z matrix.
#' #In this instance, the list contains two terms: that corresponding to the random intercept, and that to the random slope.
#' #We can also see the corresponding block of Z^T, the (transposed) random effects model matrix.
#' Ztlist<-getME(test.model,"Ztlist")
#' Ztlist
#' 
#' #we call the test using the fitted model, and indicate the random effects to be tested
#' #here we say elements in positions 2 to 2 of Ztlist are to be tested
#' #(so just element 2, the random slope, in this case)
#' 
#' lmmtest(lmmfit=test.model, Ztlist.start = 2, Ztlist.end =  2, simulations = 10**6) #should give FALSE
#' 
#' #the last two arguments are optional. 
#' #Simulations describes how many simulations we use to estimate the cdf of the weighted chisquares.
#' #Because for y, the slope effect is strong, we expect output <0.05 (H0 is rejected).
#' 
#' 
#' #Now we fit the model to y_, with weak slope effect.
#' test.model.alt<-lmer(y_~(1|f)+(0+z1|f),data=test.df)
#' summary(test.model.alt)
#' #The weak random effect is just strong enough to be detected.
#' 
#' lmmtest(test.model.alt,2,2)
#' #Because the effect is weak, we expect the output (p-value) to be >0.05.
#' 
#' #NOTE: the function will return NaN if the estimated variance of a random effect to be tested is 0.
#' #In this case H0 and H1 are the same, so it makes no sense to run the test anyway.
#'
#' 
lmmtest<-function(lmmfit, Ztlist.start, Ztlist.end, simulations=10**6){
  
  
  #Specify the positions the ranef corresponds to in the model matrix
  Ztlist<-getME(lmmfit,"Ztlist")
  pos.start=1
  pos.end=1
  if (Ztlist.start>1){
    for (index in seq(1,Ztlist.start-1) ){
      pos.start=pos.start+Ztlist[[index]]@Dim[[1]]
    }}
  for (index in seq(1,Ztlist.end) ){
    pos.end=pos.end+Ztlist[[index]]@Dim[[1]]
  }
  pos.end=pos.end-1 #we want it to be inclusive
  
  
  #Get fixed and random model matrices  
  X<-model.matrix(lmmfit,'fixed')
  Z<-model.matrix(lmmfit,'random')
  n<-nrow(X)
  M<-ncol(X)
  p<-ncol(Z)
  pm=pos.end-pos.start+1
  if (n!=nrow(Z)){
    print('Error with matrix dimensions!')
    browser()
  }
  
  #reorder Z such that the random effect is at the end.
  Z.new<-swap.column(Z,pos.start,pos.end)#fixME: I now have serious doubts about this being the right part of Z
  
  #get sigma
  sigma<-sigma(lmmfit)#should be REML fit
  
  #Get psi
  Lambdat<-getME(lmmfit,"Lambdat")
  psi<-crossprod(Lambdat) #this is the psi from the paper
  
  #swap block of psi for new ranef order
  psi<-swap.block(block.matrix=psi,pos.start = pos.start, pos.end = pos.end)
  
  #if the estimated sd of the ranef is zero, the test doesn't really make sense;
  #we return NaN as p-value.
  if (norm(psi[(p-pm+1):p,(p-pm+1):p])==0){
    return(NaN)
  }
  
  #Create psi under H0 
  zero.matrix <- Matrix(0, nrow = pm, ncol = pm, sparse = TRUE)
  psi_<-psi
  psi_[(p-pm+1):p,(p-pm+1):p]<-zero.matrix
  
  #create psi.inv and psi_.inv
  psi.inv<-solve(psi)
  psi_.inv<-pseudo_inverse(psi_)
  
  #Create B and B_
  B<-square_root(psi.inv)
  #B_<-square_root(psi_.inv) #not needed
  
  #create covariance matrix of b under H1
  cov.b<-psi*(sigma^2)
  #under H0
  cov.b_<-psi_*(sigma^2)
  
  
  
  #Assume from here on B and psi etc are correct
  
  
  #create Xtilde, augmented model matrix
  toprows<-cbind(X,Z.new)
  zero_matrix <- matrix(0, nrow = p, ncol = M)
  bottomrows<-cbind(zero_matrix,B)
  Xtilde<-rbind(toprows,bottomrows)
  
  #create P
  P<-solve(t(Xtilde)%*%Xtilde)%*%t(Xtilde)[,1:n]
  
  #create marginal covariance matrix under H0
  V_<-Z.new%*%cov.b_%*%t(Z.new) # This is Zpsi_Zt from Wood
  I<-Diagonal(n) 
  V_<-sigma^2*I+V_
  
  Pm=P[(p+M-pm+1):(p+M),]
  
  #Calculate covariance of last ranef under H0
  
  cov.bmhat<-Pm%*%V_%*%t(Pm)
  
  
  
  C<-base::chol(cov.bmhat)
  #fixME: It appears that this matrix can be rank deficient. Looking at how it is created,
  #it makes sense that P could be rank deficient, hence C being it as well.
  #If the true theoretical covariance of bmhat is rank deficient,
  #we need to cholesky factorize this semi positive def matrix anyway?.
  
  #QR decompose Xtilde, extract Rtilde 
  QR<-base::qr(Xtilde)
  Q<-qr.Q(QR)
  R<-qr.R(QR)
  Rtilde<-R[(p+M-pm+1):(p+M),(p+M-pm+1):(p+M)]
  
  #Compute the eigenvalues
  RtildeCt<-Rtilde%*%t(C)
  RtildeCt.crossprod<-crossprod(RtildeCt)
  
  e<-eigen(RtildeCt.crossprod,symmetric = TRUE, only.values = TRUE)
  eigenvalues<-e$values
  
  #Compute ||f1||^2, the test statistic
  y<-model.frame(lmmfit)[,1]
  zero.vector<-rep(0,p)
  y.0<-c(y,zero.vector)
  if (length(y.0)!=dim(Q)[1]){
    print("Error forming y.0!")
    browser()}
  
  Qty.0<-t(Q)%*%y.0
  f1<-Qty.0[(p+M-pm+1):(p+M)]
  test.statistic<-t(f1)%*%f1
  
  #perform the test
  probability<-sum.chi.squares.distr(list.weights = eigenvalues, quantile = test.statistic, simulations=simulations)
  return(1-probability) 
  #issue raised
  #I return the p-value. sum.chi.square.distr has not been updated, it can be exactly 0.
}  

# glmmtestPQL-----------------------------------------------------------------------------
#This test doesn't really follow the LMM case, instead it applies a PQL-like argument to f
#(which in this case differs from the f in Wood's book). It may work better.

#' A simple test for random effects in regression models, GLMM case.
#' @description An implementation of the test proposed by Simon Wood in the Biometrika paper "A simple test for random effects in regression models". Link: https://doi.org/10.1093/biomet/ast038
#'
#' @param glmmfit A fitted glmer() object
#' @param Ztlist.start The first index (inclusive) in Ztlist corresponding to the random effect to be tested.
#' @param Ztlist.end The last index (inclusive) in Ztlist corresponding to the random effect to be tested. For example, giving Ztlist.start=1 and Ztlist.end=2 tells lmmtest that the random effect to be tested consists of elements 1 and 2 in Ztlist. Ztlist is given by getME(lmmfit, "Ztlist"), see examples.
#' @param simulations The number of simulations used to obtain the distribution function of the test statistic. Default is 10^6.
#'
#' @return The p-value of the test. Numeric, can be exactly 0. As simulations determine the test statistic distribution repeatedly every time, results may vary slightly with each execution.
#' @export
#' @import lme4 Matrix designr
#' @rdname glmmtestPQL
#' @examples
#' #'  # For GLMMs, the test works the same as for LMMs.
#' set.seed(123)
#' 
#' #Create hierarchical dataframe
#' create_hierarchical_df<-function(levels1,levels2,obs_per_level){
#'   
#'   g1<-as.factor(1:levels1)
#'   g2<-as.factor(1:levels2)
#'   n_obs=obs_per_level #ie 20 students per class
#'   g1<-rep(g1,each=length(g2)*n_obs)
#'   g2<-rep(rep(g2,each=n_obs),times=levels1)
#'   g1<-as.factor(g1)
#'   g2<-as.factor(g2)
#'   g1.g2<-g1:g2
#'   test.df<-data.frame(g1=g1,g2=g2,g1.g2=g1.g2)
#'   return(test.df)
#'   }
#' 
#' test.df<-create_hierarchical_df(10,10,20) #function simulates a data with 10 levels of factor g1, and another 10 levels of factor g2, g2 nested within g1.
#' 
#' #Simulate the corresponding hierarchical Poisson GLMM data, with weak level 2 effect.
#' intercept<-2
#' predictor<-simLMM(~(1|g1)+(1|g1.g2),Fixef = c(intercept),VC_sd=list(c(5),c(0.1),c(0)),data=test.df)
#' y<-rpois(n=length(predictor),lambda=exp(predictor))
#' GLMM_model<-glmer(y~1+(1|g1)+(1|g1.g2),verbose=0,family = "poisson",data=test.df)
#' summary(GLMM_model)
#' 
#' 
#' getME(GLMM_model,"Ztlist")
#' #the first element of Ztlist contains the level 2 random intercept associated with the interaction factor.
#' #This is what we want to test.
#' 
#' glmmtestPQL(GLMM_model,1,1) #returns p value, or NaN if its a singular fit.

glmmtestPQL<-function(glmmfit, Ztlist.start, Ztlist.end, simulations=10**6){
  
  #Specify the positions the ranef corresponds to in the model matrix
  Ztlist<-getME(glmmfit,"Ztlist")
  pos.start=1
  pos.end=1
  if (Ztlist.start>1){
    for (index in seq(1,Ztlist.start-1) ){
      pos.start=pos.start+Ztlist[[index]]@Dim[[1]]
    }}
  for (index in seq(1,Ztlist.end) ){
    pos.end=pos.end+Ztlist[[index]]@Dim[[1]]
  }
  pos.end=pos.end-1 #we want it to be inclusive
  
  #Get fixed and random model matrices  
  X<-model.matrix(glmmfit,'fixed')
  Z<-model.matrix(glmmfit,'random')
  n<-nrow(X)
  M<-ncol(X)
  p<-ncol(Z)
  pm=pos.end-pos.start+1
  
  #reorder Z such that the random effect is at the end.
  Z.new<-swap.column(Z,pos.start,pos.end)
  
  #get sigma and W
  sigma<-sigma(glmmfit)#should be REML fit
  phi<-sigma^2
  #W_weights<-weights(glmmfit,type="working") # I think this is wrong
  sqrt_W_weights<-glmmfit@resp$sqrtWrkWt() 
  W_weights<-sqrt_W_weights^2
  W<-Diagonal(n=length(W_weights),x=W_weights)
  W.inv<-solve(W)
  sqrtW<-Diagonal(n=length(W_weights),x=sqrt_W_weights)
  
  
  #Create psi
  Lambdat<-getME(glmmfit,"Lambdat")
  psi<-crossprod(Lambdat)
  
  #swap block for new order
  psi<-swap.block(block.matrix=psi,pos.start = pos.start, pos.end = pos.end)
  
  #if the estimated sd of the ranef is zero, the test doesn't really make sense;
  #we return NaN as p-value.
  if (norm(psi[(p-pm+1):p,(p-pm+1):p])==0){
    return(NaN)
  }
  
  #Create psi_ under H0 
  
  #just replace last block by zeros
  
  zero.matrix <- Matrix(0, nrow = pm, ncol = pm, sparse = TRUE)
  psi_<-psi
  psi_[(p-pm+1):p,(p-pm+1):p]<-zero.matrix
  
  #create psi.inv and psi_.inv
  psi.inv<-solve(psi)
  psi_.inv<-pseudo_inverse(psi_)
  
  #Create B and B_
  B<-square_root(psi.inv)
  #B_<-square_root(psi_.inv) #not needed
  
  #create covariance matrix of b under H1
  cov.b<-psi*(sigma^2)
  #under H0
  cov.b_<-psi_*(sigma^2)
  
  #create Xtilde, augmented model matrix
  sqrtWX<-sqrtW%*%X
  sqrtWZ.new<-sqrtW%*%Z.new #fixME: this gives a weird warning
  toprows<-cbind(sqrtWX,sqrtWZ.new)
  zero_matrix <- base::matrix(0, nrow = p, ncol = M)
  bottomrows<-cbind(zero_matrix,B)
  Xtilde<-rbind(toprows,bottomrows)
  
  #QR decompose Xtilde, extract Rtilde 
  QR<-base::qr(Xtilde)
  Q<-qr.Q(QR) #Remark: not square/complete
  R<-qr.R(QR)
  Rtilde<-R[(p+M-pm+1):(p+M),(p+M-pm+1):(p+M)]
  
  #get pseudodata z
  #etahat<-glmmfit@resp$eta
  #e<-resid(glmmfit,type="working")
  #z<-e+etahat
  z<-glmmfit@resp$wrkResp()
  
  #get f and f1, test statistic
  zero.vector<-as.vector(rep(0,p))
  sqrtWz.0<-c(as.vector(sqrtW%*%z),zero.vector)
  if (length(sqrtWz.0)!=dim(Q)[1]){
    print("Error forming sqrtWz.0!")
    browser()}
  
  f<-t(Q)%*%sqrtWz.0
  f1<-f[(p+M-pm+1):(p+M)]
  test.statistic<-crossprod(f1)
  
  #get covariance of f under H0
  zero.matrix1<-matrix(0,nrow=n,ncol=M)
  zero.matrix2<-matrix(0,nrow=n,ncol=p)
  zero.matrix3<-matrix(0,nrow=p,ncol=M)
  toprows<-cbind(zero.matrix1,zero.matrix2)
  bottomrows<-cbind(zero.matrix3,B)
  
  B.completed<-rbind(toprows,bottomrows)
  
  psi_.completed<-matrix(0,nrow=M+p,ncol=M+p)
  psi_.completed[(M+1):(M+p),(M+1):(M+p)]<-as.matrix(psi_)
  #easier way to create these matrices, dunno why I bothered binding before
  #for some reason I have to convert psi_ to matrix, but p shouldn't be too big?
  #fixME it might become a problem
  
  
  zero.matrix1<-matrix(0,nrow = n,ncol=p)
  zero.matrix2<-matrix(0,nrow = p,ncol = p)
  toprows<-cbind(phi*Diagonal(n=n),zero.matrix1)
  bottomrows<-cbind(t(zero.matrix1),zero.matrix2)
  phiId_completed<-rbind(toprows,bottomrows)
  
  termwithR<-R-t(Q)%*%B.completed
  
  f.cov_<-t(Q)%*%phiId_completed%*%Q+phi*termwithR%*%psi_.completed%*%t(termwithR)
  f1.cov_<-f.cov_[(p+M-pm+1):(p+M),(p+M-pm+1):(p+M)]
  
  e<-eigen(f1.cov_,symmetric = TRUE, only.values = TRUE)
  eigenvalues<-e$values
  
  
  
  #get covariance of f1, under PQL normal assumption for f
  
  #removeME: these 4 could turn out irrelevant
  #psi<-crossprod(B)
  #psi_<-crossprod(B_)
  #psim<-psi[(p-pm+1):(p),(p-pm+1):(p)]
  #psim_<-psi_[(p-pm+1):(p),(p-pm+1):(p)]
  
  #f1.cov_<-phi*Diagonal(n=pm) #Under H1, a term of the type RpsiR^T gets added
  
  #get eigenvalues
  
  #removeME
  #e<-eigen(f1.cov_,symmetric = TRUE, only.values = TRUE)
  #eigenvalues<-e$values
  
  #eigenvalues<-rep(1,pm)*phi
  
  #compute test
  
  probability<-sum.chi.squares.distr(list.weights = eigenvalues, quantile = test.statistic, simulations=simulations)
  return(1-probability)
  
  
}


#Tutorial--------------------
#Tutorial: how to use the function.




# # For GLMMs, the test works the same.
# set.seed(123)
# 
# #Create hierarchical dataframe
# create_hierarchical_df<-function(levels1,levels2,obs_per_level){
#   
#   g1<-as.factor(1:levels1)
#   g2<-as.factor(1:levels2)
#   n_obs=obs_per_level #ie 20 students per class
#   g1<-rep(g1,each=length(g2)*n_obs)
#   g2<-rep(rep(g2,each=n_obs),times=levels1)
#   g1<-as.factor(g1)
#   g2<-as.factor(g2)
#   g1.g2<-g1:g2
#   test.df<-data.frame(g1=g1,g2=g2,g1.g2=g1.g2)
#   return(test.df)
#   #wasn't a return(test.df) missing before?
# }
# 
# test.df<-create_hierarchical_df(10,10,20) #function simulates a data with 10 levels of factor g1, and another 10 levels of factor g2, g2 nested within g1.
# 
# #Simulate the corresponding hierarchical Poisson GLMM data, with weak level 2 effect.
# intercept<-2
# predictor<-simLMM(~(1|g1)+(1|g1.g2),Fixef = c(intercept),VC_sd=list(c(5),c(0.1),c(0)),data=test.df)
# y<-rpois(n=length(predictor),lambda=exp(predictor))
# GLMM_model<-glmer(y~1+(1|g1)+(1|g1.g2),verbose=0,family = "poisson",data=test.df)
# summary(GLMM_model)
# 
# 
# getME(GLMM_model,"Ztlist")
# #the first element of Ztlist contains the level 2 random intercept associated with the interaction factor.
# #This is what we want to test.
# 
# glmmtestPQL(GLMM_model,1,1) #returns p value, or NaN if its a singular fit.




