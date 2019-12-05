# Title of paper:
 Alternative Models and Randomization Techniques for Bayesian Outcome Adaptive Randomization with Binary Outcomes in the ARREST Trial

# Purpose of the .R files
These programs facilitate implementing the simulation study performed in the aforementioned paper as well as generating a randomization schedule using a chosen probability model and randomization method. 

# Brief Description of the .R files
# Functions.R:
The first function ("ibb.post.prob") computes the posterior probability that the treatment arm is better than the control arm, i.e, Prob(pi_{Treatment} > pi_{Control} | Data), using the independent beta-binomial model.

The second function ("tlr.sampler") implements the logistic regression model with independent Student-t(degrees of freedom, location, scale) priors on each model parameter via a Polya-Gamma Gibbs sampler.

The third function ("tlr.post.prob") computes the posterior probability that the treatment arm is better than the control arm, i.e, Prob(pi_{Treatment} > pi_{Control} | Data), using the logistic regression model.

The fourth function ("rand.schedule") uses the posterior probability resulting from the chosen probability model to generate a randomization schedule for the next sequential group of patients using one of the three randomization methods: the weighted coin design, the mass-weighted urn design, or the modified permuted block design. This function outputs the treatment assignments for the next group of patients.

# Simulations.R
The first function ("simulate.trial") simulates the Bayesian group sequential outcome-adaptive trials that were simulated in the aforementioned paper.

The second function ("simulation.analysis") generates the desired simulations using the "simulate.trial" function and then computes the following operating characteristics: type I error rate, power, average sample size under the alternative hypothesis, average treatment arm sample size difference under the alternative hypothesis, proportion of trials in which more subjects were assigned to the inferior treatment, and root mean square error of the employed randomization method.
Users can specify the desired null response rates, sequential group size, probability model, randomization method, and the posterior probability stopping boundary, among other parameters.

The file concludes with an example of how to generate one of the simulations in the paper using these functions.

# Example.R
This file illustrates an example of how to use the functions in "Functions.R" to generate a randomization schedule for the next sequential group of patients according to the chosen trial design.
