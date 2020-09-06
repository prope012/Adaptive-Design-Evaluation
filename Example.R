# The code below provides an example of how to use the functions in Functions.R to produce a 
# randomization schedule
# The chosen probability model is the logistic regression model with a student-t prior intercept location
# of log(0.12/0.88)
# The chosen randomization method is the modified permuted block design
# accrued study data (example): 28 on treatment arm with 10 responders, 17 on control arm with 4 responders
# sequential group size of 15
# m is the number of sequential groups

setwd("C:/Users/jenni/OneDrive/Documents/Research/Bayesian Outcome Adaptive Trial Designs")
source('Functions.R')

#compute the posterior probability based on accrued study data
tlr.post.prob.results = tlr.post.prob(n=c(28,17), y=c(10,4),  m0=c(log(0.12/0.88),0))

#determine randomization schedule (for the first group, set post.prob = 0.5 instead of tlr.post.prob.results)
group.rand.schedule = rand.schedule(post.prob = tlr.post.prob.results, rand.type=c("Block"))

#1 indicates assignment to the treatment arm, 0 to the control
print(group.rand.schedule)
