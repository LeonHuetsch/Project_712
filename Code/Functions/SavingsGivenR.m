function expectAssetHoldings = SavingsGivenR(nGridAsset,minGridAsset,...
    maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt,...
    rho,r,alpha,A,depreciation,sigma,mValueGuess) 

% Given interest rate r the function calculates Savings (exptected asset
% holdings) as well as Capital and returns the absolute difference

% Capital given r

% Savings
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,...
    logShockAverage,truncOpt);


[~,mPolicyAsset,~,~] = ...
    Infinite_horizon_iteration(rho,r,alpha,A,depreciation,sigma,vGridAsset,...
    vGridShock,mTransitionShock,mValueGuess);


[~,expectAssetHoldings] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset);

%average = ((alpha*A/(r+depreciation))^(1/(1-alpha))+expectAssetHoldings)/2;
%difference = ((alpha*A/(r+depreciation))^(1/(1-alpha))-expectAssetHoldings);


end