function mixture = EMIterate(mixture, dataset)

[mixture, llnew] = EM_steps(mixture, dataset);

numSPnts = mixture.numSPnts;
D = mixture.D;
Lc = 1 + D + 0.5* D * (D+1);
epsilon = 0.01 * Lc * log(numSPnts*D);

itrNum = 1;
while true
   itrNum = itrNum + 1;
   llold = llnew;
   [mixture, llnew] = EM_steps(mixture, dataset);
   if(~mod(itrNum,10))
       fprintf('Iteration number %d ...\n llnew:%d \n',itrNum,llnew);
   end
   if (llnew - llold) <= epsilon
       break;
   end
end

mixture.loglikelilood = llnew;
mixture.epsilon = epsilon;
mixture.Lc = Lc;
