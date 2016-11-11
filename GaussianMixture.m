function [mixture] = GaussianMixture(dataset, K, ConditionNumber)
% dataset: a 1 x I cell array of scatterplots
% K: number of clusters
% mixture: will be the result of our GMM. It's a stuct with fields:
% I: number of scatterplots
% D: dim of point, for us D = 2
% K: number of clusters
% cluster(k).mu: a 1 x D array of mean 
% cluster(k).R: a D x D sigam matrix
% cluster(k).pi(i): for scatterplot i, pi is a 1 x K array of mixting coef.
% cluster(k).gamma(i): for scatterplot i, ...
%                      gamma is the znk matrix, Nd x K, ... 
%                      Nd is number of points in scatterplot i

%ConditionNumber = 1e5;
sizeDataset = size(dataset,2);
mtr = initMixture(sizeDataset, K, ConditionNumber);
% plotGMM(mtr);
mixture = EMIterate(mtr, dataset);


