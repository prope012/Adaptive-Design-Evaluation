setwd("C:/Users/jenni/Documents/T32 Research/Bayesian Outcome Adaptive Trial Designs")
source('Functions.R')

# Function to simulate Bayesian Group Sequential Outcome Adaptive Randomization trial

# response.probs = the desired response probability in each arm
# ns - analysis sample sizes
# ntrials = number of simulated trials
# max.ar - maximum randomization probability to either arm; default is 0.75
# rand.type - adaptive randomization type
#   "Coin" is biased coin flipping
#   "Urn" is a mass-weighted urn design that provides control over allocation imbalance via the max.deviation parameter
#   "Block" is a modified permuted block design
# max.deviation - maximum allowable deviation away from the targeted randomization ratio under Urn randomization
#model - probability model
#   "tlr" - logistic regression with independent t prior distributions on regression coefficients
#   "ibb" - independent beta-binomial distributions
#ibb prior specification
# pi.star - prior mean of Beta distribution, i.e. a/(a+b); default is 0.5
# pess - prior effective samplse size of Beta distribution, i.e. a+b; default is 2
#tlr prior specification
# df - prior degrees of freedom; default is 7
# m0 - prior locations; default is 0
# s0 - prior scales; default is 2.5
#tlr mcmc sampler specification
# nsamps = the number of Polya-Gamma Gibbs Sampler iterations; default is 2000
# warmup = the number of iterations in the burn-in period (tossed out) for the Polya-Gamma Gibbs Sampler; default is 100

simulate.trial <- function(response.probs, ns, max.ar=0.75, rand.type=c("Coin","Urn","Block"), max.deviation=3, model=c("tlr","ibb"),
                           pi.star=0.5, pess=2, df=c(7,7), m0=c(0,0), s0=c(2.5,2.5), nsamps=2000, warmup=100){
  
  #Initialize Data
  n = c(0,0); y = c(0,0); rand.prob = 0.5
  
  #Storage object for posterior probabilities
  stats = matrix(NA,nrow=length(ns),ncol=5)
  colnames(stats) = c("PP","nE","yE","nC","yC")
  rownames(stats) = 1:length(ns)
  
  #Generate Data under Block AR scheme and true response probabilities
  for(group in 1:nrow(stats)){
    #Number of new enrollees during current group
    n.new = c(0,ns)[group+1]-c(0,ns)[group]
    
    #Generate assignments and outcomes for current interval
    if(rand.type == "Coin"){
      nE.new = rbinom(1,n.new,rand.prob)
      yE.new = rbinom(1,nE.new,response.probs[1])
      nC.new = n.new - nE.new
      yC.new = rbinom(1,nC.new,response.probs[2])
    }
    if(rand.type == "Urn"){
      z = NULL; rand.prob.temp = rand.prob
      for(i in 1:n.new){
        z[i] = rbinom(1,1,min(max(rand.prob.temp,0),1))
        rand.prob.temp = rand.prob.temp+(rand.prob-1*(z[i]==1))/max.deviation
      }
      nE.new = sum(z); nC.new = n.new - nE.new
      yE.new = rbinom(1,nE.new,response.probs[1])
      yC.new = rbinom(1,nC.new,response.probs[2])
    }
    if(rand.type == "Block"){
      u = rbinom(1,1,n.new*rand.prob - floor(n.new*rand.prob))
      nE.new = u*ceiling(n.new*rand.prob)+(1-u)*floor(n.new*rand.prob)
      nC.new = n.new-nE.new
      yE.new = rbinom(1,nE.new,response.probs[1])
      yC.new = rbinom(1,nC.new,response.probs[2])
    }
    #Update dataset
    n = n+c(nE.new,nC.new)
    y = y+c(yE.new,yC.new)
    
    #Update posterior probability that CCL is superior to ICU and randomization probability for next group
    if(model=="ibb") post.prob = ibb.post.prob(n=n,y=y,pi.star=pi.star,pess=pess)
    if(model=="tlr") post.prob = rowMeans(tlr.sampler(n=n,y=y,X=rbind(c(1,0.5),c(1,-0.5)),df=df,m0=m0,s0=s0,nsamps=nsamps,warmup=warmup)>0)[2]
    
    rand.prob = min(max.ar,max(1-max.ar,post.prob))
    
    #Store Relevant Statistics
    stats[group,] = c(post.prob,n[1],y[1],n[2],y[2])  
  }
  
  return(stats)
}

# ns defines the sequential group size for each trial as well as the maximum sample size 
# A set.seed of 1994 was used for all simulations
# ntrials = the desired number of simulations; we used 10000, but the below code uses 100
# null.response.probs = the null response rates evaluated for each simulated scenario
# The code below uses a group size of 15 and a maximum sample size of 150; however, we also evaluated
# the below statistics at group sizes of 30 and 50
# The code below is ran using the IBB model with prior mean 0.50 and Weighted Coin randomizion method; 
# however, we also evaluated the below statistics using the 3 other probability models and 2 other randomization methods
#   The other models: ibb with prior mean 0.12, logistic regression with prior intercept locations of 
#   0 and log(0.12/0.88)
#   The other randomization methods: modified permuted block design and the mass weighted urn design
# When implementing the logistic regression model we recommend using a parallel computation package to 
# increase compuation speed (unparallelized code below)
# PP.bound = posterior probability stopping boundary; each unique model and group size combination
# require different stopping boundaries - these were chosen to ensure that power and type I error (12%
# null response rates) were approximately 0.90 and 0.05, respectively
# nE/nC = number of subjects assigned to the treatment arm/control arm

ntrials = 10000
ns = seq(15,150,15) #group size of 15 for a maximum sample size of 150
null.response.probs = c(0.04,0.08,0.12,0.16,0.20,0.30,0.40,0.50,0.75,0.95)

#Function to run simulations and compute desired statistics
simulation.analysis = function(null.response.probs, ns, max.ar=0.75, rand.type=c("Coin","Urn","Block"), max.deviation=3, model=c("tlr","ibb"),
                               pi.star=0.5, pess=2, df=c(7,7), m0=c(0,0), s0=c(2.5,2.5), nsamps=2000, warmup=100, PP.bound){
  
  # run simulations at the desired null response rates and the hypothesized alternative response rates
  Results = list()
  for (i in 1:length(null.response.probs)){
    Results[[i]] = lapply(1:ntrials,function(trial) simulate.trial(response.probs = c(null.response.probs[i],null.response.probs[i]), ns=ns, max.ar=max.ar, rand.type=rand.type, 
                         max.deviation=max.deviation, model=model, pi.star=pi.star, pess=pess, df=df, m0=m0, s0=s0, nsamps=nsamps, warmup=warmup))
  }
  
  # compute type I error using inputed posterior probability stopping boundary
  bounds= c(rep(PP.bound,length(ns)-1),PP.bound)  
  Type.I.Error = data.frame(rep(0,length(Results)))
  for (i in 1:length(Results)){Type.I.Error[i,] = mean(sapply(Results[[i]],function(x) max(x[,1] > bounds | x[,1] < 1-bounds) == 1))}
  Type.I.Error = cbind(null.response.probs,Type.I.Error)
  colnames(Type.I.Error) = c("Null Response Probabilities","Type I Error")
  
  
  # compute power using inputed posterior probability stopping boundary and hypothesized alternative
  # response rates of 0.37 for treatment and 0.12 for control
  Alt.Scenario = lapply(1:ntrials,function(trial) simulate.trial(response.probs = c(0.37,0.12), ns=ns, max.ar=max.ar, rand.type=rand.type, max.deviation=max.deviation, model=model, 
                                                                 pi.star=pi.star, pess=pess, df=df, m0=m0, s0=s0, nsamps=nsamps, warmup=warmup))
  Power = data.frame(mean(sapply(Alt.Scenario,function(x) max(x[,1] > bounds | x[,1] < 1-bounds) == 1)))
  colnames(Power) = "Power"
  
  # compute average overall sample size (n), average treatment arm sample size difference (nE-nC), and 
  # proportion of nE < nC, i.e, proportion of trials with more individuals assigned to the inferior
  # treatment, under the alternative hypothesis after applying posterior probability stopping boundaries
  
  alt.foo.15 = sapply(Alt.Scenario,function(x){
    stop.at = which.max(x[,1] > c(bounds[-length(ns)],0)) #determine when trial exceeds posterior probability stopping boundary for the first time
    n = sum(x[stop.at,c(2,4)]) #total number of people in trial before stopping, i.e., nE+nC
    nE = x[stop.at,2] #total number of people on treatment arm before stopping
    nC = x[stop.at,4] #total number of people on control arm before stopping
    nE.minus.nC = nE - nC #treatment arm sample size difference
    nE.less.nC = nE<nC
    
    return(c(n,nE.minus.nC,nE.less.nC)) #returns desired parameters: n, nE-nC, and if nE-nC < 0
  })
  Alt.Results <- t(data.frame(rowMeans(alt.foo.15))); #returns average n, nE-nC, and proportion of trials in which nE-nC<0
  colnames(Alt.Results) <- c("Average Sample Size", "Average nE-nC","mean(nE-nC<0)")
  rownames(Alt.Results) = NULL
  
  # compute RMSE for the employed randomization method, probability model, and group size at the 
  # hypothesized null response rate of 12% and the alternative hypothesis response rates
 
  # function to generate Squared Error
  SE = function(x){
    stop.at = which.max(x[,1] > c(bounds[-length(ns)],0)) #determine when trial exceeded PP stopping boundary
    if (stop.at == 1){rand.probs = 0.5} #a randomization probability of 0.5 was used to generate treatment assignments for first sequential group
    else {
      rand.probs = c(0.5,x[1:(stop.at-1),1]) #determine randomization probabilities used to generate treatment assignments in each sequential group
          for (b in 1:length(rand.probs)){ #max and min randomization probabilities are 0.75 and 0.25
               if (rand.probs[b] < 0.25){rand.probs[b] = 0.25}
               if (rand.probs[b] > 0.75){rand.probs[b] = 0.75}}}
    expected.nE = sum(rand.probs*(150/length(ns))) #determine expected number of subjects assigned to treatment group
    actual.nE = x[stop.at,2] #determine actual number of subjects assinged to treatment group
    nE.Error = expected.nE - actual.nE #calculate difference between expected and actual number of treatment assignments
    SE = (nE.Error)^2 #compute Squared Error
    return(SE)
    }
  
  #compute RMSE of interest
  RMSE.of.interest <- list(Results[[3]],Alt.Scenario) #scenarios generated from hypothesized null (12%) and alternative (37%, 12%) response rates
  RMSE.results = data.frame(rep(0,2)) #create object to store RMSE results
  for(i in 1:2) {RMSE.results[i,] = sqrt(mean(sapply(RMSE.of.interest[[i]], function(x) SE(x))))}

  rownames(RMSE.results) = c("null", "alternative")
  colnames(RMSE.results) = "RMSE"
  
  return(c(list(Type.I.Error, Power, Alt.Results, RMSE.results)))
  
}

# Use function to compute above statistics for IBB model with prior mean 0.50, group size 15, and
# weighted coin randomization
set.seed(1994)
Example = simulation.analysis(null.response.probs = null.response.probs, rand.type = "Coin", model = "tlr", pi.star = 0.50, PP.bound = 0.9872, ns=ns)