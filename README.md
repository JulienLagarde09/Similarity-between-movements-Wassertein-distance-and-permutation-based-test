# Similarity-between-movements-Wasserstein-distance-and-permutation-based-test
Test on two samples distributions based on Wassertein distance and permutations

function [observdistance, threshold] = testwasspermutation(x,y,alpha,permutations)

Aim: Testing for difference between two processes using cdf ("fonction de
répartition in french", càd integral(histo(x))).
The test is non parametric: Permutation.
The idea is simple: compute the distance between cdf of actual (observed) data,
then to get a null H0: exchange randomly data between samples and compute distance.
NB the distance recalls KS testing on max distance between cdf. Here use
Wasserstein stuff, or EMD, = the sum of differences between cdf.
To estimate the distance uses this implementation: emd_dist.m,
to estimate Wasserstein distance using cumsum(hist(x)). binning maybe not
optimal.
For the permutations uses: randperm([x;y])
can be done with exact total N permutations (ouch... O(n?))
distance is positive, hence one we do a one sided test
threshold = (100 - alpha) percentiles of
distribution of N distances between randomized samples from actual data,
The null H0 = random exchangeability.
Permutation for null H: solid if iid (?; check).
distance > 0, if observ > threshold reject the null H0
INPUTS: data in x and y, one column vector each, alpha = sig level (0.05 or 0.01 for ex),
permutations = for example 100 randomizations, use depends on length of the data
OUTPUTS: the observed distance, and the threshold under the null
a sample of references:
in R, twosamples package (uses various distances)
https://cran.r-project.org/web/packages/twosamples/readme/README.html
https://arxiv.org/abs/2007.01360
Ramdas, A., Trillos, N. G., & Cuturi, M. (2017). On wasserstein two-sample testing and
related families of nonparametric tests. Entropy, 19(2), 47.
http://www.gatsby.ucl.ac.uk/~gretton/mmd/mmd.htm
permutation (randomization): see Sir R Fisher, 'course :)
NB: tests exist using L2 Wasserstein distance
Julien Lagarde, Uni Montpellier, Euromov dec 2020
