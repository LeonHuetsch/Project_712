function [mStationaryDist,expectAssetHoldings] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset)


% Calculates the stationary distribution. It first calculates the 
% transition matrix for the states with the transition matrix of 
% the shocks (income process) and policy function for asset holdings
 
nGridAsset = length(vGridAsset);
nStates = nGridShock*nGridAsset;

% Transition matrix (a1y1,a1y2,..., a2y1,a1y2,...) for measure Phi
mTransitionState=zeros(nStates);
for ac=1:nGridAsset
    for yc=1:nGridShock
        ind=find(vGridAsset==mPolicyAsset(ac,yc));
        mTransitionState(yc+(ac-1)*nGridShock,(ind-1)*nGridShock...
            +1:(ind-1)*nGridShock+nGridShock) = mTransitionShock(yc,:);
    end
end

% First eigenvector has eigenvalue 1, iteration is faster though
%[mEigenvector,vEigenvalue] = eig(mTransitionState');
%vEigenvector1 = mEigenvector(:,1);
%vStationaryDist = vEigenvector1/sum(vEigenvector1);

vStationaryDist = ones(nStates,1)/(nStates);
maxit=1e04;
it=1;
tol=1e-06;
diff=100;

while it<maxit && diff>tol
    vHelp = mTransitionState' * vStationaryDist;
    
    diff = max(abs(vHelp-vStationaryDist));
    vStationaryDist = vHelp;
    it=it+1;
end

mStationaryDist = reshape(vStationaryDist,[nGridShock,nGridAsset])';
expectAssetHoldings = sum(sum(mStationaryDist.*mPolicyAsset));

end