




library(lme4)
library(Matrix)
library(designr)

set.seed(123)

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

# Create a test block matrix to test swap.block, swap.column on
M1 <- matrix(1:16, nrow = 4, ncol = 4)
M2 <- matrix(1:9 , nrow=3)
M3 <- matrix(1:4, nrow=2)
test.matrix<- bdiag(M1,M2,M3)
test.matrix.swap<-as.matrix(swap.block(test.matrix,5,7))
as.matrix(swap.column(test.matrix,5 ,7))



create_hierarchical_df<-function(levels1,levels2,obs_per_level){
  
  g1<-as.factor(1:levels1)
  g2<-as.factor(1:levels2)
  n_obs=obs_per_level #ie 20 students per class
  g1<-rep(g1,each=length(g2)*n_obs)
  g2<-rep(rep(g2,each=n_obs),times=levels1)
  g1<-as.factor(g1)
  g2<-as.factor(g2)
  g1.g2<-g1:g2
  test.df<-data.frame(g1=g1,g2=g2,g1.g2=g1.g2)
  return(test.df)
  #wasn't a return(test.df) missing before?
}

# glmmtestPQL--------------
#This test doesn't really follow the LMM case, instead it applies a PQL-like argument to f
#(which in this case differs from the f in Wood's book). It may work better.

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
  
  #determine whether the null model still contains ranefs
  no.ranef.null=FALSE
  if(Ztlist.end-Ztlist.start+1==length(Ztlist)){
    no.ranef.null=TRUE
  }
  
  #Get fixed and random model matrices  
  X<-model.matrix(glmmfit,'fixed')
  Z<-model.matrix(glmmfit,'random')
  n<-nrow(X)
  M<-ncol(X)
  p<-ncol(Z)
  pm=pos.end-pos.start+1
  
  #reorder Z such that the random effect is at the end.
  Z.new<-swap.column(Z,pos.start,pos.end)#fixME: I now have serious doubts about this being the right part of Z Update: no this should be dandy
  
  #get sigma and W
  sigma<-sigma(glmmfit)#should be REML fit
  phi<-sigma^2
  #W_weights<-weights(glmmfit,type="working") # I think this is wrong
  sqrt_W_weights<-glmmfit@resp$sqrtWrkWt() 
  W_weights<-sqrt_W_weights^2
  W<-Diagonal(n=length(W_weights),x=W_weights)
  W.inv<-solve(W)
  sqrtW<-Diagonal(n=length(W_weights),x=sqrt_W_weights)
  
  
  #Create block diagonal B, under H1
  B<-getME(glmmfit,"Lambdat")
  #suffices to swap block of B to get cholesky factor for new ranef order
  B<-swap.block(block.matrix=B,pos.start = pos.start, pos.end = pos.end)
  Bt<-t(B)
  #if the estimated sd of the ranef is zero, the test doesn't really make sense;
  #we return NaN as p-value.
  if (norm(B[(p-pm+1):p,(p-pm+1):p])==0){
    return(NaN)
  }
  
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
  bottomrows<-cbind(zero.matrix3,B_)
  
  B_completed<-rbind(toprows,bottomrows)
  
  psi.completed<-matrix(0,nrow=M+p,ncol=M+p)
  psi_<-crossprod(B_)
  psi.completed[(M+1):(M+p),(M+1):(M+p)]<-as.matrix(psi_)
  #easier way to create these matrices, dunno why I bothered binding before
  #for some reason I have to convert psi_ to matrix, but p shouldn't be too big?
  #fixME it might become a problem
  
  
  zero.matrix1<-matrix(0,nrow = n,ncol=p)
  zero.matrix2<-matrix(0,nrow = p,ncol = p)
  toprows<-cbind(phi*Diagonal(n=n),zero.matrix1)
  bottomrows<-cbind(t(zero.matrix1),zero.matrix2)
  phiId_completed<-rbind(toprows,bottomrows)
  
  termwithR<-R-t(Q)%*%B_completed
  
  f.cov_<-t(Q)%*%phiId_completed%*%Q+phi*termwithR%*%psi.completed%*%termwithR
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



#glmmtest variacnce of z-------------------
glmmtest<-function(glmmfit,Ztlist.start,Ztlist.end,simulations=10**6,H_0=TRUE){
  
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
  
  #determine whether the null model still contains ranefs
  no.ranef.null=FALSE
  if(Ztlist.end-Ztlist.start+1==length(Ztlist)){
    no.ranef.null=TRUE
  }
  
  #Get fixed and random model matrices  
  X<-model.matrix(glmmfit,'fixed')
  Z<-model.matrix(glmmfit,'random')
  n<-nrow(X)
  M<-ncol(X)
  p<-ncol(Z)
  pm=pos.end-pos.start+1
  
  #reorder Z such that the random effect is at the end.
  Z.new<-swap.column(Z,pos.start,pos.end)#fixME: I now have serious doubts about this being the right part of Z Update: no this should be dandy
  
  #get sigma and W
  sigma<-sigma(glmmfit)#should be REML fit
  phi<-sigma^2
  #W_weights<-weights(glmmfit,type="working") # I think this is wrong
  sqrt_W_weights<-glmmfit@resp$sqrtWrkWt() 
  W_weights<-sqrt_W_weights^2
  W<-Diagonal(n=length(W_weights),x=W_weights)
  W.inv<-solve(W)
  sqrtW<-Diagonal(n=length(W_weights),x=sqrt_W_weights)
  
  #Create block diagonal B, under H1
  B<-getME(glmmfit,"Lambdat")
  #suffices to swap block of B to get cholesky factor for new ranef order
  B<-swap.block(block.matrix=B,pos.start = pos.start, pos.end = pos.end)
  Bt<-t(B)
  #if the estimated sd of the ranef is zero, the test doesn't really make sense;
  #we return NaN as p-value.
  if (norm(B[(p-pm+1):p,(p-pm+1):p])==0){
    return(NaN)
  }
  
  #Create B under H0 
  
  #just replace last block of B by zeros
  
  zero.matrix <- Matrix(0, nrow = pm, ncol = pm, sparse = TRUE)
  B_<-B
  B_[(p-pm+1):p,(p-pm+1):p]<-zero.matrix
  B_t<-t(B_)
  
  #create Xtilde, augmented model matrix
  sqrtWX<-sqrtW%*%X
  sqrtWZ.new<-sqrtW%*%Z.new #fixME: this gives a weird warning
  toprows<-cbind(sqrtWX,sqrtWZ.new)
  zero_matrix <- base::matrix(0, nrow = p, ncol = M)
  bottomrows<-cbind(zero_matrix,B)
  Xtilde<-rbind(toprows,bottomrows)
  
  #create P and Pm
  P<-solve(t(Xtilde)%*%Xtilde)%*%t(Xtilde)[,1:n]
  Pm=P[(p+M-pm+1):(p+M),]
  
  
  #QR decompose Xtilde, extract Rtilde 
  QR<-base::qr(Xtilde)
  Q1<-qr.Q(QR) #Remark: not square/complete
  R<-qr.R(QR)
  Rtilde<-R[(p+M-pm+1):(p+M),(p+M-pm+1):(p+M)]
  
  #get pseudodata z
  z<-glmmfit@resp$wrkResp() #fixMe: Is this correct? it was wrong before
  #get f and f1, test statistic
  zero.vector<-as.vector(rep(0,p))
  sqrtWz.0<-c(as.vector(sqrtW%*%z),zero.vector)
  if (length(sqrtWz.0)!=dim(Q1)[1]){
    print("Error forming sqrtWz.0!")
    browser()}
  
  f<-t(Q1)%*%sqrtWz.0 #actually not true, this is f and r together
  f1<-f[(p+M-pm+1):(p+M)]
  test.statistic<-crossprod(f1)
  
  #get covariance of f1, under PQL normal assumption for f
  
  psi<-crossprod(B)
  psi_<-crossprod(B_)
  b.cov_<-phi*psi_
  b.cov<-phi*psi
  psim<-psi[(p-pm+1):(p),(p-pm+1):(p)]
  psim_<-psi_[(p-pm+1):(p),(p-pm+1):(p)]#well this should be zero, dunno why thats here#removeME
  
  #sqrtWz.cov_<-phi*Diagonal(n=n)+sqrtWZ.new%*%b.cov_%*%t(sqrtWZ.new)
  sqrtWz.cov_<-phi*Diagonal(n=n)+phi*tcrossprod(sqrtWZ.new%*%B_t) #might be better numerically
  #fixME: If I put b.cov instead of b.cov_, I do in fact get reasonable p values.
  if (H_0==FALSE){
    sqrtWz.cov_<-phi*Diagonal(n=n)+phi*tcrossprod(sqrtWZ.new%*%Bt) #is sqrtWz.cov, without the _, really #removeME
  }
  #However, under H0 the test seems unreasonably powerful (often p value of 0, but sometimes like 0.6).. weird
  
  bmhat.cov_<-Pm%*%sqrtWz.cov_%*%t(Pm)
  C<-base::chol(bmhat.cov_)
  
  #get eigenvalues
  RtildeCt<-Rtilde%*%t(C)
  RtildeCt.crossprod<-crossprod(RtildeCt)
  
  e<-eigen(RtildeCt.crossprod,symmetric = TRUE, only.values = TRUE)
  eigenvalues<-e$values
  
  #compute test
  probability<-sum.chi.squares.distr(list.weights = eigenvalues, quantile = test.statistic, simulations=simulations)
  return(1-probability)
  
  
}

#Create hierarchical dataframe
#test.df<-create_hierarchical_df(20,10,5)
test.df<-create_hierarchical_df(10,10,20) 

#simulate Poisson GLMM
intercept<-4
predictor<-simLMM(~(1|g1)+(1|g1.g2),Fixef = c(intercept),VC_sd=list(c(5),c(0.02),c(0)),data=test.df)
y<-rpois(n=length(predictor),lambda=exp(predictor))
GLMM_model<-glmer(y~0+(1|g1)+(1|g1.g2),offset = rep(intercept,length(y)),verbose=0,family = "poisson",data=test.df)
summary(GLMM_model)



#Do the test
debug(glmmtest)
glmmtest(GLMM_model,1,1)
undebug(glmmtest)

#Checking weights W
GLMM_model@resp$sqrtWrkWt()^2-weights(GLMM_model,"working")
#oh okay so this is the same thing, perhaps I should switch back to the working weights then.



#Using the function------------
#Could the problem have to do with offsets? Do it again, this time let glmer estimate the intercept, no offsets

intercept<-2
predictor<-simLMM(~(1|g1)+(1|g1.g2),Fixef = c(intercept),VC_sd=list(c(5),c(0.1),c(0)),data=test.df)

y<-rpois(n=length(predictor),lambda=exp(predictor))
GLMM_model<-glmer(y~1+(1|g1)+(1|g1.g2),verbose=0,family = "poisson",data=test.df)
summary(GLMM_model)

debug(glmmtest)
glmmtest(GLMM_model,1,1,H_0=TRUE) #the real test
glmmtest(GLMM_model,1,1,H_0=FALSE) #comparing with H1
undebug(glmmtest)

getME(GLMM_model,"Ztlist")

#Simulate the p-values under H1 to see if the function is reasonable under H1, and the problem lies with the modified covariance of hatb under H0.
simulations<-100
pvalues<-rep(NaN,simulations)
pvalues_<-rep(NaN,simulations)
for (i in 1:simulations){
  intercept<-4
  predictor<-simLMM(~(1|g1)+(1|g1.g2),Fixef = c(intercept),VC_sd=list(c(5),c(0.05),c(0)),verbose=FALSE,data=test.df)
  y<-rpois(n=length(predictor),lambda=exp(predictor))
  GLMM_model<-glmer(y~1+(1|g1)+(1|g1.g2),verbose=0,family = "poisson",data=test.df)
  
  result<-tryCatch({
    # Code to be executed
    #pvalues[i]<-glmmtest(GLMM_model,1,1,simulations=5*10**4,H_0=FALSE)
    pvalues_[i]<-glmmtest(GLMM_model,1,1,simulations=10**3,H_0=TRUE)
  }, error = function(err) {
    # Code to handle the error
    pvalues[i] <- NaN
  })
  
}
plot(sort(pvalues))
abline(a=0,b=0.01)
#p values are reasonable, slightly larger than expected


plot(sort(pvalues_))
abline(a=0,b=0.01)
sort(pvalues_)
#this is suspicious


#Shall we try a power simulation anyway?
simulations<-50
denominator<-400
p_values<-rep(NaN,simulations)
sd_vector<-1:simulations/denominator
for (i in 1:simulations){
  sd<-i/denominator 
  predictor<-simLMM(~(1|g1)+(1|g1.g2),Fixef = c(intercept),VC_sd=list(c(5),c(sd),c(0)),verbose=FALSE,data=test.df)
  y<-rpois(n=length(predictor),lambda=exp(predictor))
  tryCatch({
  GLMM_model<-glmer(y~1+(1|g1)+(1|g1.g2),verbose=0,family = "poisson",data=test.df)
  p_values[i]<-glmmtest(GLMM_model,1,1,simulations=10**4)
  },error=function(err){  }
  )
}#end for
plot(sd_vector,p_values)
p_values


hierarchical_simulation<-function(VC_sd=list(c(5),c(0),c(0)),nr_tests=10**2,significance_level=0.05){
  
  H0_count<-0
  H1_count<-0
  p_value_vector<-c()
  for (testnr in 1:nr_tests){
    print(paste("testnr ",testnr))
    predictor<-simLMM(~(1|g1)+(1|g1.g2),Fixef = c(intercept),VC_sd=list(c(5),c(sd),c(0)),verbose=FALSE,data=test.df)
    y<-rpois(n=length(predictor),lambda=exp(predictor))
    p_value<-NaN
    tryCatch({
      GLMM_model<-glmer(y~1+(1|g1)+(1|g1.g2),verbose=0,family = "poisson",data=test.df)
      p_value<-glmmtest(GLMM_model,1,1,simulations=10**4)
    }, error=function(err){p_value<-NaN}) #the error thing doesnt seem to work
    p_value_vector<-append(p_value_vector,p_value)
    if(is.nan(p_value)){
      H0_count<-H0_count+1
    } else if (p_value>=0.05){
      H0_count<-H0_count+1
    } else {H1_count<-H1_count+1}
    
  }#end for
  result=list(H0_count,H1_count,p_value_vector)
  return(result)
}
#Power for some different values of sd of the tested ranef

sd_list<-0.01*(1:15)
power_list<-function(sd_list,nr_tests=3){
  result=rep(0,length(sd_list))
  for (counter in 1:length(sd_list)){
    output<-hierarchical_simulation(VC_sd = list(5,sd_list[[counter]],0),nr_tests = nr_tests)
    H1_count<-output[[2]]
    p_value_vector<-output[[3]]
    print(p_value_vector)
    print(H1_count)
    result[counter]<-mean(p_value_vector)#modif, now this is p value
  }
  return(result)
}
power<-power_list(sd_list = sd_list)
power
plot(sd_list,100*power,main="Power Curve",xlab="Deviation of Level 2 Random Intercept",ylab="Power (%)")

#Power Analysis glmmPQL----
hierarchical_simulationPQL<-function(VC_sd=list(c(5),c(0),c(0)),nr_tests=10**2,significance_level=0.05){
  
  H0_count<-0
  H1_count<-0
  p_value_vector<-c()
  for (testnr in 1:nr_tests){
    print(paste("testnr ",testnr))
    predictor<-simLMM(~(1|g1)+(1|g1.g2),Fixef = c(intercept),VC_sd=VC_sd,verbose=FALSE,data=test.df)
    y<-rpois(n=length(predictor),lambda=exp(predictor))
    p_value<-NaN
    tryCatch({
      GLMM_model<-glmer(y~1+(1|g1)+(1|g1.g2),verbose=0,family = "poisson",data=test.df)
      p_value<-glmmtestPQL(GLMM_model,1,1,simulations=10**4)
    }, error=function(err){
      p_value<-NaN
      print(err)}) #the error thing doesnt seem to work
    #Why did I need a try again? Does the model not converging generate an error?
    #Oh I believe for strange fits/nonconvergence the cholesky factorization fails sometimes.
    
    p_value_vector<-append(p_value_vector,p_value)
    if(is.nan(p_value)){
      H0_count<-H0_count+1
    } else if (p_value>=0.05){
      H0_count<-H0_count+1
    } else {H1_count<-H1_count+1}
    
  }#end for
  result=list(H0_count,H1_count,p_value_vector)
  return(result)
}

result<-hierarchical_simulationPQL(VC_sd=list(c(5),c(0.05),c(0)))
mean(result[[3]],na.rm=TRUE)
result[[2]]/100 # power 
#results:
#sd 0.2 : approx 0.06 p value
#sd 0.05 approx 0.33

#more systematically:
nr_points=15 #nr of points on the graph
powerPQL<-rep(0,nr_points)
p_valuesPQL<-rep(0,nr_points)
for (i in 1:nr_points){
  sd<-i*0.02
  result<-hierarchical_simulationPQL(VC_sd=list(c(5),c(sd),c(0)))
  powerPQL[i]<-result[[2]]/100
  p_valuesPQL[i]<-mean(result[[3]],na.rm=TRUE)
}
plot((1:nr_points)*0.02,p_valuesPQL,xlab="SD of tested Random Effect",ylab = "Average p value",main = "Performance of glmmtestPQL")
plot((1:nr_points)*0.02,powerPQL,xlab="SD of tested Random Effect",ylab="Power",main = "Performance of glmmtestPQL") #fantastic results!!!

#Significance Test

result<-hierarchical_simulationPQL(nr_tests=400)

result[[1]]/400#significance level 
result[[3]]# p-value

#apparently, no tests rejected H0, and the p values are abnormally concentrated around 0.5. Worrying.
