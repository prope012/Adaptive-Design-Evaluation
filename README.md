# Title:
 Alternative Models and Randomization Techniques for Bayesian Outcome Adaptive Randomization with Binary Outcomes in the ARREST Trial

# Purpose of the .R files
These programs facilitate implementing the simulation study performed in the aforementioned paper as well as generating a randomization schedule using a chosen probability model and randomization method. 

# Brief Description of the .R files
# Functions.R:
The first function ("ibb.post.prob") computes the posterior probability that the treatment arm is better than the control arm, i.e, Prob(pi_{Treatment} > pi_{Control} | Data), using the independent beta-binomial model.

The second function ("tlr.sampler") implements the logistic regression model with independent Student-t(degrees of freedom, location, scale) priors on each model parameter via a Polya-Gamma Gibbs sampler.

The third function ("tlr.post.prob") computes the posterior probability that the treatment arm is better than the control arm, i.e, Prob(pi_{Treatment} > pi_{Control} | Data), using the logistic regression model.

The fourth function ("rand.schedule") uses the posterior probability resulting from the chosen probability mdoel to generate a randomization schedule for the next sequential group of patients using one of the three randomization methods: the weighted coin design, the mass-weighted urn design, or the modified permuted block design. This function outputs the treatment assignments for the next group of patients.

# Simulations.R

# Example.R

