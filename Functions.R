# The functions below compute the Posterior Probability that the Treatment Arm is better than the 
# Control Arm, i.e., Prob(pi_{Treatment} > pi_{Control} | Data), using one of the following 3 
# probability models: Independent Beta-Binominal (ibb) model, Logit Transformation Model, and Logistic 
# Regression Model

# Below is the function to compute the Posterior Probability that the Treatment Arm is better than the 
# Control Arm with the ibb Model, that is yE/yC ~ Binomial(nE/nC,pi_E/piC) and pi_E/pi_C ~ Beta(pi0*pess,(1-pi0)*pess)
# where yE/yC = number of responders on the treatment arm/number of responders of the control arm,
# pi_E/pi_C = response rate of the treatment arm/response rate of the control arm
# pi.star = prior mean of Beta distribution, i.e. a/(a+b); default is 0.50
# pess = prior effective samplse size of Beta distribution, i.e. a+b; default is 2
# n = a vector of size 2 where the first and second numbers indicate the number of subjects assigned to
# the treatment and control arms, respectively
# y = a vector of size 2 where the first and second numbers indicate the number of successes on the
# treatment and control arms, respectively

ibb.post.prob <- function(n,y,pi.star=0.5,pess=2){
  post.prob <- unlist(integrate(function(x) pbeta(x,y[1]+pi.star*pess,(n[1]-y[1])+(1-pi.star)*pess,lower.tail=FALSE)*dbeta(x,y[2]+pi.star*pess,(n[2]-y[2])+(1-pi.star)*pess),lower=0,upper=1))$value
  return(post.prob)
}

# The first function below implements the logit transformation and logisitc regression models with 
# independent t(degrees of freedom, location, scale) priors on each model parameter via a Polya-Gamma 
# Gibbs sampler
# To implement the logit transformation model, input a design matrix of X = rbind(c(1,0),c(0,1)) into
# the tlr.sampler function (the default)
# To implement the logistic regression model, input a design matrix of X = rbind(c(1,0.5),c(1,-0.5))
# The second function computes the Posterior Probability that the Treatment Arm is better than the 
# Control Arm using the tlr.sampler function; ensure that you specify which model (the logit transformation
# model or logistic regression model) the posterior probability is to be calculated for
# n = vector of size 2; first and second numbers indicate the number of subjects assigned to the 
# treatment and control arms, respectively (default is 1)
# y = vector of size 2; first and second numbers indicate the number of successes on the treatment 
# and control arms, respectively
# df = the prior degrees of freedom for the regression coefficients; default is 7 
# m0 = the prior locations for the regression coefficients; default is 0
# s0 = the prior scales for the regression coefficients; default is 2.5
# nsamps = the number of Polya-Gamma Gibbs Sampler iterations; default is 2000
# warmup = the number of iterations in the burn-in period (tossed out) for the Polya-Gamma Gibbs 
# Sampler; default is 100
# tlogit = logit transformation model with independent Student-t priors on the model coefficients
# tlr = logistic regression model with independent Student-t priors on the regression coefficients

library(BayesLogit)
tlr.sampler = function(n=NULL,y=c(1,0),X=rbind(c(1,0),c(0,1)),df=rep(7,ncol(X)),m0=rep(0,ncol(X)),s0=rep(2.5,ncol(X)),nsamps=2000,warmup=100){
  
  #All trials are of size 1
  if(is.null(n)) n = rep(1,length(y))
  
  # Pseudo-outcome vector
  kappa = cbind(y-0.5*n)
  
  # Storage object
  beta.samps = matrix(NA,nrow=ncol(X),ncol=nsamps)
  
  # Initial values
  beta = rep(0,ncol(X))
  
  # Gibbs Sampler 
  for(samp in 1:(nsamps+warmup)){
    
    # Sample beta
    omega = rpg.devroye(ncol(X),n,as.vector(X%*%cbind(beta))) # generate auxillary PG parameters, i.e. omegas
    tau = rgamma(ncol(X),shape=(df+1)/2,rate=(df*s0^2+(beta-m0)^2)/2) # generate auxillary Gamma parameters, i.e. taus
    V = chol2inv(chol(t(X)%*%(X*omega) + diag(tau))) # Posterior variance matrix
    m = V%*%(t(X)%*%kappa+diag(tau)%*%cbind(m0)) # Posterior mean vector
    beta = as.vector(m + t(chol(V))%*%rnorm(ncol(X))) # generate regression parameters, i.e. beta
    
    # Save posterior sample 
    if(samp>warmup){
      beta.samps[,samp-warmup] = beta
    } 
  }
  
  # Return posterior samples
  return(beta.samps) 
}

#Compute posterior probability based on the tlr.sampler function; input study data for n and y
tlr.logit.post.prob <- function(model = c("tlr","tlogit"),n,y,df=c(7,7),m0=c(0,0),s0=c(2.5,2.5),nsamps=2000,warmup=100){
  
  #if chosen probability model is the logistic regression model
  if(model=="tlr") post.prob = rowMeans(tlr.sampler(n=n,y=y,X=rbind(c(1,0.5),c(1,-0.5)),df=df,m0=m0,s0=s0,nsamps=nsamps,warmup=warmup)>0)[2]
  
  #if chosen probability model is the logit transformation model
  if(model=="tlogit") {post.prob = mean((tlr.sampler(n=n,y=y,X=rbind(c(1,0),c(0,1)),df=df,m0=m0,s0=s0,nsamps=nsamps,warmup=warmup)[1,]) >
                      (tlr.sampler(n=n,y=y,X=rbind(c(1,0),c(0,1)),df=df,m0=m0,s0=s0,nsamps=nsamps,warmup=warmup)[2,]))}
  return(post.prob)
  }

# The below function uses the posterior probability computed above via the chosen probability model to
# generate a randomization schedule for the next sequential group of patients
# Three different randomization methods are available as an input: the weighted coin design
# (Coin), the Mass-Weighted Urn Design (Urn) and the Modified Permuted Block Design (Block)
# group.size = size of sequential group; by default, the function below uses this number as the block 
# size for the MPBD randomization method; default is 15
# max.ar = the maximum randomization probability allowed to either arm; default is 0.75
# maximum deviation = maximum allowable deviation away from the targeted randomization ratio under MWUD
# randomization; default is 3
# nE.new = number of new subjects to be assigned to the treatment arm
# nC.new = number of new subjects to be assigned to the control arm
# post.prob = the posterior probability used to generate the randomization schedule; corresponds to
# the chosen probability model (ibb, tlogit, or tlr)
# Outputs the treatment assignments for the next group of patients

rand.schedule = function(post.prob, group.size = 15, max.ar=0.75, rand.type=c("Coin","Urn","Block"), max.deviation=3) {

  #determine the randomization probability for the treatment group
  rand.prob = min(max.ar,max(1-max.ar,post.prob))

  #generate number of subjects assigned to each treatment via the chosen randomization method
  if(rand.type == "Coin"){
  nE.new = rbinom(1,group.size,rand.prob)
  nC.new = group.size - nE.new
  }
  if(rand.type == "Urn"){
  z = NULL; rand.prob.temp = rand.prob
  for(i in 1:group.size){
    z[i] = rbinom(1,1,min(max(rand.prob.temp,0),1))
    rand.prob.temp = rand.prob.temp+(rand.prob-1*(z[i]==1))/max.deviation
  }
  nE.new = sum(z)
  nC.new = group.size - nE.new
  }
  if(rand.type == "Block"){
  u = rbinom(1,1,group.size*rand.prob - floor(group.size*rand.prob))
  nE.new = u*ceiling(group.size*rand.prob)+(1-u)*floor(group.size*rand.prob)
  nC.new = group.size-nE.new
  }

  #generate group randomization schedule
  schedule = sample(c(rep(1,nE.new),rep(0,nC.new)))
  return(schedule) 
}
