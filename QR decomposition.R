# Understanding the QR decomposition
library(Matrix)

#creating a LMM, we will do QR decomp on the augmented model matrix
number.factors<-100 
f<-1:number.factors
f<-as.factor(f)
z1<-rnorm(5*number.factors)
intercept<-3
test.df<-data.frame(f=f,
                    z1=z1)

y<-simLMM(formula=~(1|f)+(0+z1|f),Fixef = c(intercept),VC_sd=list(c(5),c(5),c(2)),data=test.df)

y_<-simLMM(formula=~(1|f)+(0+z1|f),Fixef = c(intercept),VC_sd=list(c(5),c(1),c(2)),data=test.df)

test.df$y<-y
test.df$y_<-y_
test.model<-lmer(y~(1|f)+(0+z1|f),data=test.df)
test.model_<-lmer(y~(1|f),data=test.df)

test.model.alt<-lmer(y_~(1|f)+(0+z1|f),data=test.df)
test.model.alt_<-lmer(y_~(1|f),data=test.df)

#create augmented model matrix
X<-model.matrix(test.model.alt,"fixed")
Z<-model.matrix(test.model.alt,"random")
toprows<-cbind(X,Z)
B<-getME(test.model.alt,"Lambdat")
p<-dim(Z)[2]
M<-dim(X)[2]
n<-dim(X)[1]
zero.matrix<-Matrix(0,ncol=M,nrow=p,sparse = TRUE)
bottomrows<-cbind(zero.matrix,B)
Xtilde<-rbind(toprows,bottomrows)

#do naive qr decomposition from matrix package
QR<-qr(Xtilde,backPermute=TRUE)
Q<-qr.Q(QR)
R<-qr.R(QR)
R.alt<-qrR(QR)
R.mat<-as.matrix(R)
R.alt.mat<-as.matrix(R.alt)

Q%*%R #gives wrong matrix

Q%*%R.alt #gives correct matrix

#Analyzing QR objects
str(QR)
dim(QR@R)#this must give the full R matrix
dim(R)#this gives the restricted square block
norm(QR@R[1:201,]-R)#indeed, this is 0.

#It turns out qr returns an annoyingly permutated result
#PAP_ = QR, with qr.Q outputting (P_)^{-1}Q  

q<-QR@q
p<-QR@p
#p and q encode the permutations

#Supposedly PA=A[p+1,]
#then we should have PA=qr.Q(..)R

norm(Xtilde[,q+1]-Q%*%R)#e-14, so basically 0

#To Recap: qr computes PAP_=QR, and it returns qr.Q=P^{-1}Q, qr.R=R, and qrR=R(P_)^{-1}.

#QR@p gives us a vector representing the permutation p, that is to say 
#A[p+1,]=PA
#while QR@q represents P_, that is to say
#A[,q+1]=AP_

#Proof:

norm(Q[p+1,]%*%R-(Xtilde[p+1,])[,q+1]) #basically 0, proving that QR=PAP_


#This function should invert a permutation (starting at 1, not 0)
invert.permutation<-function (vector){
  n<-length(vector)
  result<-rep(0,n)
  for (i in 1:n){
    j<-vector[i]
    result[j]<-i
  }
  return(result)
}

#let's see if the same happens for dense matrices, can we circumvent the problem by converting to dense matrices?
Xtilde.mat<-as.matrix(Xtilde)
QR.converted<-qr(Xtilde.mat)
Q.converted<-qr.Q(QR.converted)
R.converted<-qr.R(QR.converted)
norm(Xtilde.mat-Q.converted%*%R.converted) #oh wow it seems to work, but obviously numerically its not good.
Q.converted%*%R.converted
Xtilde

#Using the base function qr also works, though seems to have the same numerical problems
QR.base<-base::qr(Xtilde)
Q.base<-base::qr.Q(QR.base)
R.base<-base::qr.R(QR.base)
norm(R.base-R)
norm(Q.base%*%R.base-Xtilde)
Q.base%*%R.base-Xtilde
Xtilde


#base::qr even works for rank-deficient matrices, where the decomposition is not unique.
M=matrix(c(1,-1,-1,1),nrow=2,ncol=2)
qr<-base::qr(M)
Q<-qr.Q(qr)
R<-qr.R(qr)
Q
R
Q%*%R
