function emdv1v2 = emd_dist(v1,v2)
% use : emdv1v2 = emd_dist(v1,v2)
% principle: distance = sum of differences between cdf of stochastic processes
% cdf = integral of pdf (cumsum of hist)
% Wasserstein or earth mover distance. Theory in C Villani 2003
% C. Villani (2003), Topics in Optimal Transportation, American Mathematical Society
% Muskulus, M., Houweling, S., Verduyn-Lunel, S., & Daffertshofer, A. (2009).
% Functional similarities and distance properties. Journal of neuroscience methods, 183(1), 31-41.
% Code from: in Piotr et al. paper supp material
% Slowinski, P.& Tsaneva-Atanasova, K. (2016). Dynamic similarity promotes interpersonal 
% coordination in joint action. Journal of The Royal Society Interface, 13(116), 20151093.
min_v1 = min(v1); min_v2 = min(v2);% get the range of data (min & max)
max_v1 = max(v1); max_v2 = max(v2);
min_all = min(min_v1,min_v2);
max_all = max(max_v1,max_v2);

bins = linspace(min_all,max_all,101); % support Z with 101 bins
binwidth = bins(2)- bins(1); % widths of the bins, i.e. dz
maxemd = abs(max_all-min_all); % maximal EMD

h1 = hist(v1,bins); % h1 and h2 are velocity profiles
h2 = hist(v2,bins); % v1 and v2 are velocity time series
l1 = numel(v1); % for normalisation
l2 = numel(v2); % for normalisation
emdv1v2 = sum(abs(cumsum(h1/l1)-cumsum(h2/l2)))*binwidth/maxemd;