function [mixture,likelihood] = EM_steps(mixture,dataset)

I = mixture.I;
K = mixture.K;
D = mixture.D;
numPnts = 0;
numSPnts = 0;
likelihood = 0;
prePik = [];
preMu = [];

mu = double(zeros(K,D));
Nik = double(zeros(I,K));
for i = 1:I
    pixels = double(dataset{i}.pixels);
    Nd = size(pixels,1);
    numPnts = numPnts + Nd;
    w = double(dataset{i}.w);
    h = double(dataset{i}.h);
    salVal = double(dataset{i}.salValue);
    N = sum(salVal);
    numSPnts = numSPnts + N;
    pixels = 100*([pixels(:,1)/w,pixels(:,2)/h]);
    znk = zeros(Nd,K);
    
    % First compute posteriors (EStep) and find likelihood
    for k =1:K
        Y1 = pixels - ones(Nd,1) * mixture.cluster(k).mu;
        Y2 = -0.5 * Y1*mixture.cluster(k).invR;
        znk(:,k) = dot(Y1,Y2,2) + mixture.cluster(k).const;
        pik(k) = mixture.cluster(k).pi(i);
    end
    prePik(i,:) = pik;
    llmax = max(znk,[],2);
    znk = exp(znk - llmax * ones(1,K));
    znk = znk .* (ones(Nd,1) * pik);
    ss = sum(znk,2);
    znk = znk ./ (ss*ones(1,K));
    likelihood = likelihood + sum(salVal .* (log(ss) + llmax));

    % Now compute unnormailized-mu and pi
    for k=1:K
       mu(k,:) = mu(k,:) + (salVal .* znk(:,k))' * pixels; 
       Nik(i,k) = sum(salVal .* znk(:,k));
       % compute pi_ik
       mixture.cluster(k).pi(i) = Nik(i,k)/N;
    end   
end
for k = 1:K
    Nk(k) = sum(Nik(:,k));     
    mu(k,:) = mu(k,:)/Nk(k);
    preMu(k,:) = mixture.cluster(k).mu;
    mixture.cluster(k).mu = mu(k,:);
end

% Now compute sigma
for k = 1:K
    R{k} = double(zeros(D,D));
end
for i = 1:I
    pixels = double(dataset{i}.pixels);
    Nd = size(pixels,1);
    w = double(dataset{i}.w);
    h = double(dataset{i}.h);
    salVal = double(dataset{i}.salValue);
    pixels = 100*([pixels(:,1)/w,pixels(:,2)/h]);
    znk = zeros(Nd,K);
    
    % First compute posteriors 
    for k =1:K
        Y1 = pixels - ones(Nd,1) * preMu(k,:);
        Y2 = -0.5 * Y1*mixture.cluster(k).invR;
        znk(:,k) = dot(Y1,Y2,2) + mixture.cluster(k).const;
    end
    pik = prePik(i,:);
    llmax = max(znk,[],2);
    znk = exp(znk - llmax * ones(1,K));
    znk = znk .* (ones(Nd,1) * pik);
    ss = sum(znk,2);
    znk = znk ./ (ss*ones(1,K));
    
    for k = 1:K
        for r = 1:D
            for s = r:D
                Ri(r,s) = ...
                    ((pixels(:,r)-mixture.cluster(k).mu(r))' * ...
                    ((pixels(:,s)-mixture.cluster(k).mu(s)) .* ...
                    (znk(:,k) .* salVal)));
                if r ~= s
                    Ri(s,r) = Ri(r,s);
                end           
            end
        end
        R{k} = R{k} + Ri;
    end
end

for k = 1:K
    mixture.cluster(k).R = R{k}/Nk(k) + mixture.Rmin*eye(mixture.D);
    mixture.cluster(k).invR = inv(mixture.cluster(k).R);
    mixture.cluster(k).const = -(D * log(2*pi) + ...
                               log(det(mixture.cluster(k).R)))/2;
end
mixture.numPnts = numPnts;
mixture.numSPnts = numSPnts;


