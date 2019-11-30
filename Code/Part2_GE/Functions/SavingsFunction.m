function savings = SavingsFunction(vGridAsset,vMultiSteps,vGridShock,mTransitionShock,nGridShock,rrho,r,aalpha,A,depreciation,ssigma,mValFunGuess) 

% Given interest rate r the function calculates Savings (exptected asset
% holdings) as well as Capital and returns the absolute difference

% Capital given r

minAsset = vGridAsset(1);
maxAsset = vGridAsset(end);
[~,~,mPolicyAsset,~] = MultigridVFI_InfHorizon(rrho,r,ssigma,aalpha,A,depreciation,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,mValFunGuess,10);
%    VFI_InfHorizon(rrho,r,ssigma,aalpha,A,depreciation,vGridAsset,...
%    vGridShock,mTransitionShock,mValFunGuess,10);

[~,savings] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset);

%average = ((alpha*A/(r+depreciation))^(1/(1-alpha))+expectAssetHoldings)/2;
%difference = ((alpha*A/(r+depreciation))^(1/(1-alpha))-expectAssetHoldings);

end