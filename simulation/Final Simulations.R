library(lme4)
library(Matrix)
library(designr)
#library(Woodsimplevartest)

set.seed(123)


#Analysis of lmmtest-------------

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

#this function simulates data and records the p values, the times H0 wasn't rejected, the times H0 was rejected.
#for hierarchical data of exactly the type of test.df (for illustration purposes)
#Arguments: variances, number of times to run the function and record results, nr simulations to be used in lmmtest
#Output: list of length 3 containing p values, #H0 not rejected, # H0 rejected.
hierarchical_simulation_lmm<-function(VC_sd=list(c(5),c(0),c(2)),nr_tests=10**2,significance_level=0.05){
  
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

#Significance testing
set.seed(123)

result<-hierarchical_simulation(nr_tests = 400)
result[[2]]/400 #Type 1 Error rate
mean(result[[3]],na.rm=TRUE) #average p value, NaNs removed

#Simulating Power 

#we vary the standard deviation of the ranef to be tested
sd_list<-0.2*(1:15)
power_list<-function(sd_list,nr_tests=10**2){
  result=rep(0,length(sd_list))
  for (counter in 1:length(sd_list)){
    output<-hierarchical_simulation_lmm(VC_sd = list(5,sd_list[[counter]],2),nr_tests = nr_tests)
    H1_count<-output[[2]]
    result[counter]<-H1_count/nr_tests
  }
  return(result)
}
power<-power_list(sd_list = sd_list)
plot(sd_list,100*power,main="Performance of lmmtest",xlab="Deviation of Level 2 Random Intercept",ylab="Power (%)")



#Analysis of glmmtestPQL-----------------------------

#the equivalent function to hierarchical_simulation_lmm, for glmmtestPQL.
#The standard deviations to be given are now on the predictor level. Notably, the last value should always be zero.
#Default is for testing significance level, ie level 2 ranef is zero as well.
hierarchical_simulationPQL<-function(VC_sd=list(c(5),c(0),c(0)),nr_tests=10**2,significance_level=0.05){
  
  H0_count<-0
  H1_count<-0
  p_value_vector<-c()
  for (testnr in 1:nr_tests){
    print(paste("testnr ",testnr))
    predictor<-simLMM(~(1|g1)+(1|g1.g2),Fixef = c(intercept),VC_sd=VC_sd,verbose=FALSE,data=test.df)
    p_value<-"Error"
    while(p_value=="Error"){
      try({
        y<-rpois(n=length(predictor),lambda=exp(predictor))
        GLMM_model<-glmer(y~1+(1|g1)+(1|g1.g2),verbose=0,family = "poisson",data=test.df)
        p_value<-glmmtestPQL(GLMM_model,1,1,simulations=10**4)
      })
    }
    #Why did I need a try again? 
    #Rare cases of model not converging or Cholesky factorization failing.
    
    p_value_vector<-append(p_value_vector,p_value)
    if(is.nan(p_value)){
      H0_count<-H0_count+1
    } else if (p_value>=0.05){
      H0_count<-H0_count+1
    } else if(p_value<0.05){H1_count<-H1_count+1}
    
  }#end for
  result=list(H0_count,H1_count,p_value_vector)
  return(result)
}

# Significance Test

result<-hierarchical_simulationPQL(nr_tests=400)

result[[1]]/400#significance level 
result[[3]]# p-value

#apparently, no tests rejected H0, and the p values are concentrated around 0.5. Unusual.

# Power
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
