function [mStationaryDist,savings] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset)


% Calculates the stationary distribution. It first calculates the 
% transition matrix for the states with the transition matrix of 
% the shocks (income process) and policy function for asset holdings
 
nGridAsset = length(vGridAsset);
nStates = nGridShock*nGridAsset;


maxit=1e04;
it=1;
tol=1e-06;
diff=100;
mStationaryDist = ones(nGridAsset,nGridShock)/nStates ;

while it<maxit && diff>tol
    it = it+1;
    newDist = zeros(nGridAsset,nGridShock);
    for shock = 1:nGridShock
        for cap=1:nGridAsset
          newDist(mPolicyAsset(cap,shock),:) = newDist(mPolicyAsset(cap,shock),:) + mStationaryDist(cap,shock)*mTransitionShock(shock,:);
        end
    end  
    diff=max(max(abs(newDist-mStationaryDist)));
    mStationaryDist = newDist;
end

savings = sum(sum(mStationaryDist.*vGridAsset(mPolicyAsset)));
end