function expectAssetHoldings = SavingsGivenR(vGridAsset,vGridShock,...
    mTransitionShock,nGridShock,rrho,r,aalpha,A,depreciation,ssigma,mValueGuess) 

% Given interest rate r the function calculates Savings (exptected asset
% holdings) as well as Capital and returns the absolute difference

% Capital given r

[~,~,mPolicyAsset,~] = ...
    VFI_InfHorizon(rrho,r,ssigma,aalpha,A,depreciation,vGridAsset,...
    vGridShock,mTransitionShock,mValueGuess,10);


[~,expectAssetHoldings] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset);

%average = ((alpha*A/(r+depreciation))^(1/(1-alpha))+expectAssetHoldings)/2;
%difference = ((alpha*A/(r+depreciation))^(1/(1-alpha))-expectAssetHoldings);


end