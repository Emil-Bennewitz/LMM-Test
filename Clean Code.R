library(lme4)
library(Matrix)
library(designr)

set.seed(123)

### cdf of the weighted sum of chisquare random variables

#Note: total distribution is simulated "simulations" times <-increasing load with nr factor levels/ranefs 
#gets too much when hundreds of factor levels

sum.chi.squares.distr<-function(list.weights,quantile,simulations=1000000){
  sum=rep(0,simulations)
  for (weight in list.weights){
    sum=sum+weight*rchisq(simulations,1)
  }
  result<-ecdf(sum)(quantile)
  return(result)
}

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

# Create a test block matrix to test swap.block, swap.column on
M1 <- matrix(1:16, nrow = 4, ncol = 4)
M2 <- matrix(1:9 , nrow=3)
M3 <- matrix(1:4, nrow=2)
test.matrix<- bdiag(M1,M2,M3)
test.matrix.swap<-as.matrix(swap.block(test.matrix,5,7))
as.matrix(swap.column(test.matrix,5 ,7))


##### ---------------------------------------------------lmmtest------------------------------------

#fixME: weights are currently ignored for simplicity
lmmtest<-function(lmmfit, Ztlist.start, Ztlist.end, no.ranef.null=FALSE, simulations=10**6){
  
  
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
  Z.new<-swap.column(Z,pos.start,pos.end)
  
  #get sigma
  sigma<-sigma(lmmfit)#should be REML fit
  
  #Create block diagonal B, under H1
  B<-getME(lmmfit,"Lambdat")
  #suffices to swap block of B to get cholesky factor for new ranef order
  B<-swap.block(block.matrix=B,pos.start = pos.start, pos.end = pos.end)
  Bt<-t(B)
  B.matrix<-as.matrix(B)
  
  #Create B under H0 
  
  #just replace last block of B by zeros
  if (no.ranef.null==FALSE){
    zero.matrix <- Matrix(0, nrow = pm, ncol = pm, sparse = TRUE)
    B_<-B
    B_[(p-pm+1):p,(p-pm+1):p]<-zero.matrix
    B_t<-t(B_)
    
  }#end if
  
  #create covariance matrix of b, psi
  #under H1
  psi<-crossprod(B)*(sigma^2)
  #under H0
  if (no.ranef.null==FALSE){
    psi_<-crossprod(B_)*(sigma^2)}
  
  
  #create Xtilde, augmented model matrix
  toprows<-cbind(X,Z.new)
  zero_matrix <- matrix(0, nrow = p, ncol = M)
  bottomrows<-cbind(zero_matrix,B)
  Xtilde<-rbind(toprows,bottomrows)
  
  #create P
  P<-solve(t(Xtilde)%*%Xtilde)%*%t(Xtilde)[,1:n]
  
  #create marginal covariance matrix under H0
  if (no.ranef.null==FALSE){
    V_<-tcrossprod(Z.new%*%B_t) # This is Zpsi_Zt from Wood
    I<-Diagonal(n) 
    V_<-sigma^2*(I+V_)} #now its the marginal cov
  else {
    I<-Diagonal(n)
    V_<-sigma^2*I
  }
  
  Pm=P[pos.start:pos.end,]
  
  #Calculate covariance of last ranef under H0
  
  cov.bmhat<-Pm%*%V_%*%t(Pm)
  C<-chol(cov.bmhat)
  
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
  return(probability)
}  

#Tutorial: how to use the function.

#Example: LMM, one factor, 100 levels, 5 observations per factor.
#Random intercept, random slope, y=3+alpha+z*beta+epsilon, as follows:
number.factors<-100 
f<-1:number.factors
f<-as.factor(f)
z1<-rnorm(5*number.factors)
intercept<-3
test.df<-data.frame(f=f,
                    z1=z1)
#simulate responses: y has a strong random slope effect, y_ a very weak effect
y<-simLMM(formula=~(1|f)+(0+z1|f),Fixef = c(intercept),VC_sd=list(c(5),c(5),c(2)),data=test.df)
y_<-simLMM(formula=~(1|f)+(0+z1|f),Fixef = c(intercept),VC_sd=list(c(5),c(0.4),c(2)),data=test.df)

test.df$y<-y
test.df$y_<-y_

#fit models to the simulated data
test.model<-lmer(y~(1|f)+(0+z1|f),data=test.df)


#looking at this object, we can see the different components of the Z matrix.
#In this instance, the list contains two terms: that corresponding to the random intercept, and that to the random slope.
#We can also see the corresponding block of Z^T, the (transposed) random effects model matrix.
Ztlist<-getME(test.model,"Ztlist")
Ztlist

#we call the test using the fitted model, and indicate the random effects to be tested
#here we say elements in positions 2 to 2 of Ztlist are to be tested
#(so just element 2, the random slope, in this case)

lmmtest(lmmfit=test.model, Ztlist.start = 2, Ztlist.end =  2, no.ranef.null = FALSE, simulations = 10**6) #should give FALSE

#the last two arguments are optional. 
#no.ranef.null will be removed, but indicates whether any random effects are left under H0,
#because the code is slightly different in this case.
#Simulations describes how many simulations we use to estimate the cdf of the weighted chisquares.
#Because for y, the slope effect is strong, we expect output >0.95 (H0 is rejected).


#Now we fit the model to y_, with weak slope effect.
test.model.alt<-lmer(y_~(1|f)+(0+z1|f),data=test.df)
summary(test.model.alt)
#The weak random effect is just strong enough to be detected.

lmmtest(test.model.alt,2,2)
#Because the effect is weak, we expect the output to be <0.95.

#NOTE: the function will fail if the estimated variance of a random effect to be tested is 0.
#This is because the augmented model matrix, Xtilde, is rank deficient.
#The matrix P=(XtildeT*Xtilde)^{-1}XtildeT cannot be calculated.
#However, if a random effect is estimated to be 0 by lmer(), we should obviously exclude it,
#and there is no need to test for it. 
#This is why I chose y to have a strong and y_ to have a weak, but present, random slope effect.
