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
  
  #determine whether the null model still contains ranefs
  no.ranef.null=FALSE
  if(Ztlist.end-Ztlist.start+1==length(Ztlist)){
    no.ranef.null=TRUE
  }
  
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
  
  #Create block diagonal B, under H1
  B<-getME(lmmfit,"Lambdat")
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
  if (no.ranef.null==FALSE){#fixME: are these ifs even necessary?
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
    V_<-sigma^2*(I+V_)} else {
    I<-Diagonal(n)
    V_<-sigma^2*I
  }
  
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

lmmtest(lmmfit=test.model, Ztlist.start = 2, Ztlist.end =  2, simulations = 10**6) #should give FALSE

#the last two arguments are optional. 
#no.ranef.null will be removed, but indicates whether any random effects are left under H0,
#because the code is slightly different in this case.
#Simulations describes how many simulations we use to estimate the cdf of the weighted chisquares.
#Because for y, the slope effect is strong, we expect output <0.05 (H0 is rejected).


#Now we fit the model to y_, with weak slope effect.
test.model.alt<-lmer(y_~(1|f)+(0+z1|f),data=test.df)
summary(test.model.alt)
#The weak random effect is just strong enough to be detected.

lmmtest(test.model.alt,2,2)
#Because the effect is weak, we expect the output (p-value) to be >0.05.

#NOTE: the function will fail if the estimated variance of a random effect to be tested is 0.
#This is because the augmented model matrix, Xtilde, is rank deficient.
#The matrix P=(XtildeT*Xtilde)^{-1}XtildeT cannot be calculated.
#However, if a random effect is estimated to be 0 by lmer(), we should obviously exclude it,
#and there is no need to test for it. 
#This is why I chose y to have a strong and y_ to have a weak, but present, random slope effect.




#New example: Hierarchical models.

g1<-as.factor(1:10)
g2<-as.factor(1:5)
n_obs=20 #ie 20 students per class
g1<-rep(g1,each=length(g2)*n_obs)
g2<-rep(rep(g2,each=n_obs),times=10)
g1<-as.factor(g1)
g2<-as.factor(g2)
g1.g2<-g1:g2
test.df<-data.frame(g1=g1,g2=g2,g1.g2=g1.g2)

y<-simLMM(~(1|g1)+(1|g1.g2),Fixef=c(0),VC_sd=list(c(5),c(4),c(2)),data=test.df)
#Somewhat clumsy way to simulate hierarchical data, I couldn't find better;
#(1|g1)+(1|g1:g2) didn't work, so I just created the third factor g1.g2, 
#If you know how to use simLMM for hierarchical models, please help

y_<-simLMM(~(1|g1)+(1|g1.g2),Fixef=c(0),VC_sd=list(c(5),c(0.4),c(2)),data=test.df)

test.model<-lmer(y~0+(1|g1)+(1|g1:g2),data=test.df)
test.model.alt<-lmer(y_~0+(1|g1)+(1|g1:g2),data=test.df)
summary(test.model.alt)

Ztlist<-getME(test.model,"Ztlist")
Ztlist
length(Ztlist)

lmmtest(test.model,1,1)

lmmtest(test.model.alt,1,1) #for deviance 0.4, I got 0.62, for dev=1, 0.99.



#Analyzing the significance level-------------

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
}

test.df<-create_hierarchical_df(10,10,20)

#simulate the test multiple times and count how many times H0 (the given sd) is true, how many times H1.
hierarchical_simulation<-function(VC_sd=list(c(5),c(0),c(2)),nr_tests=10**2,significance_level=0.05){

  H0_count<-0
  H1_count<-0
  p_value_vector<-c()
  for (testnr in 1:nr_tests){
    print(paste("testnr ",testnr))
    y_<-simLMM(~(1|g1)+(1|g1.g2),Fixef=c(0),VC_sd=VC_sd,data=test.df,verbose=FALSE)
    test.model.alt<-lmer(y_~0+(1|g1)+(1|g1:g2),data=test.df)
    p_value<-lmmtest(test.model.alt,1,1,simulations=10**5)#10^5 sims to run faster
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
debug(hierarchical_simulation)
counts<-hierarchical_simulation(VC_sd = list(5,1,2))
undebug(hierarchical_simulation)
H0_count<-counts[[1]]
H1_count<-counts[[2]]
p_value_vector<-counts[[3]]
H0_count
H1_count
sort(p_value_vector)


#Power for some different values of sd of the tested ranef

sd_list<-0.2*(1:15)
power_list<-function(sd_list,nr_tests=10**2){
  result=rep(0,length(sd_list))
  for (counter in 1:length(sd_list)){
    output<-hierarchical_simulation(VC_sd = list(5,sd_list[[counter]],2),nr_tests = nr_tests)
    H1_count<-output[[2]]
    result[counter]<-H1_count/nr_tests
  }
  return(result)
}
power<-power_list(sd_list = sd_list)
plot(sd_list,100*power,main="Performance of lmmtest",xlab="Deviation of Level 2 Random Intercept",ylab="Power (%)")

power

#Significance Testing

result<-hierarchical_simulation(nr_tests = 400)
result[[2]]/400
mean(result[[3]],na.rm=TRUE)
result[[3]]

#playground

y_<-simLMM(~(1|g1)+(1|g1.g2),Fixef=c(0),VC_sd=list(c(5),c(0),c(2)),data=test.df)
test.model.alt<-lmer(y_~0+(1|g1)+(1|g1:g2),data=test.df)
summary(test.model.alt)
getME(test.model.alt,"Ztlist")
norm(getME(test.model.alt,"Lambdat")[1:50,1:50])#the corresponding block of B is zero.
ranef(test.model.alt)
lmmtest(test.model.alt,1,1)#NaNs now implemented


