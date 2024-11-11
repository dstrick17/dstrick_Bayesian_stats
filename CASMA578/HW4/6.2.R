# 6.2 Mixture model:
"The file glucose.dat contains the plasma glucose concentration of 532 females 
 from a study on diabetes (see Exercise 7.6)."

glucose_conentration <-  data <- read.table("path/to/your/data.txt", header = FALSE, sep = "", stringsAsFactors = FALSE)


"a) Make a histogram or kernel density estimate of the data. Describe
 how this empirical distribution deviates from the shape of a normal
 distribution"




"b) Consider the following mixture model for these data: For each study
participant there is an unobserved group membership variable Xi
which is equal to 1 or 2 with probability p and 1 − p. If Xi = 1
then Yi ∼ normal(θ1, σ1**U2), and if Xi = 2 then Yi ∼ normal(θ2, σ2**2). Let
p ∼ beta(a, b), θj ∼ normal(µ0,τ0)**2 and 1/σj ∼ gamma(ν0/2, ν0σ0**2/2
for both j = 1 and j = 2. Obtain the full conditional distributions of
(X1, . . . , Xn), p, θ1, θ2, σ1**2 and σ2**2"


"c) Setting a = b = 1, µ0 = 120, τ0**2= 200, σ0**2 = 1000 and ν0 = 10,
implement the Gibbs sampler for at least 10,000 iterations. Let
θ1**(s) = min{θ1**(s), θ2**(s)} and θ2**(s) = max{θ1**(s), θ2**(s)}. 
Compute and plot the autocorrelation functions of θ1**(s) and θ2**(s), 
as well as their effective sample sizes."


"d) For each iteration s of the Gibbs sampler, sample a value x ∼
binary(p**(s)), then sample Y˜(s) ∼ normal(θx**(s), σx**2(s). 
Plot a histogram or kernel density estimate for the empirical distribution of
Y˜**(1), . . . , Y˜**(S), and compare to the distribution in part a). 
Discuss the adequacy of this two-component mixture model for the glucose
data."