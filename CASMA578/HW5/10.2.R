# 10.2 Nesting success: Younger male sparrows may or may not nest during a
# mating season, perhaps depending on their physical characteristics. Researchers have recorded the nesting success of 43 young male sparrows
# of the same age, as well as their wingspan, and the data appear in the
# file msparrownest.dat. Let Yi be the binary indicator that sparrow i
# successfully nests, and let xi denote their wingspan. Our model for Yi
# is logit Pr(Yi = 1|α, β, xi) = α + βxi, where the logit function is given by
# logit θ = log[θ/(1 − θ)].

# a) Write out the joint sampling distribution Qn
# i=1 p(yi|α, β, xi) and simplify as much as possible.

# b) Formulate a prior probability distribution over α and β by considering the range of Pr(Y = 1|α, β, x) as x ranges over 10 to 15, the
# approximate range of the observed wingspans.

# c) Implement a Metropolis algorithm that approximates p(α, β|y, x).
# Adjust the proposal distribution to achieve a reasonable acceptance
# rate, and run the algorithm long enough so that the effective sample
# size is at least 1,000 for each parameter.

# d) Compare the posterior densities of α and β to their prior densities.

# e) Using output from the Metropolis algorithm, come up with a way to
# make a confidence band for the following function fαβ(x) of wingspan:
#   fαβ(x) = e α+βx1 + eα+βx ,
# where α and β are the parameters in your sampling model. Make a
# plot of such a band.


# Load data set
sparrow_data <- read.table("https://raw.githubusercontent.com/dstrick17/dstrick_Bayesian_stats/refs/heads/main/CASMA578/HW5/msparrownest.txt", header=TRUE)
