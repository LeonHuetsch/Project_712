tic;

addpath('VFI_Matrix')
addpath('Functions')
outpath ='Output/Output_Part2/';

r = 0.02;
rho = 0.04;
alpha = 0.66;
depreciation = 1;
%A = 1;
A = @(r) (r+depreciation)/alpha;

delta = 0.8;   % Persistence of income shock
sigmaY = 0.2;   % Variance of income shock    

nGridAsset = 100;
nGridShock = 21;

minGridAsset = 0;
maxGridAsset = 10;

logShockAverage = 0;
truncOpt = 0;


[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt);


[mValueFunction,mPolicyAsset,~,~] = ...
    Infinite_horizon_iteration(rho,r,vGridAsset,...
    vGridShock,mTransitionShock);


[mStationaryDist,expectAssetHoldings] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset);

%capStock = capitalDemand(interestRateRCE);


toc;