function mixture = initMixture(datasetSize, K, ConditionNumber)

mixture.I = datasetSize;
mixture.K = K;
mixture.D = 2;

mu = load('initMuMat');
mu = mu.mu;
R = load('initCovMat');
R = R.sigma;
mixture.Rmin = mean(diag(R))/ConditionNumber;

pik = 1/K * ones(datasetSize,1);
for k = 1:K
   cluster(k).mu = mu(k,:); 
   cluster(k).R = R + mixture.Rmin * eye(mixture.D);
   cluster(k).pi = pik;
   cluster(k).invR = inv(cluster(k).R);
   cluster(k).const = -(mixture.D * log(2*pi) + log(det(cluster(k).R)))/2;
end
mixture.cluster = cluster;

