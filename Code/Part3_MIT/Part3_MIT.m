%% Housekeeping

close all
addpath('Functions')
outpath ='Output/';




%% Parameters

bbeta = 0.96;
rrho = 1/bbeta - 1;
aalpha = 0.36;
depr = 0.08;

ssigma = 1;
ddelta = 0.8;
ssigmaY = 0.2;

llambda = 0;
ttau = 0;
kkappa = 0.5;

nGridAsset = 250;
nGridShock = 21;
vMultiSteps = [50,250];
optAccelerator = 10;  

minAsset = 0;
maxAsset = 40;

logShockAverage = 0;
truncOpt = 0;




%% Stationary equilibrium for A=1

tic
A = 1;
    
EqCondition = @(r) abs(ConditionsGE(nGridAsset,minAsset,maxAsset,vMultiSteps,nGridShock,ssigmaY,ddelta,...
    logShockAverage,truncOpt,rrho,aalpha,A,depr,ssigma,0,llambda,r,kkappa,ttau,optAccelerator));

%options = optimset('fminbnd');
%options.TolX = 1e-03;

[rSS,diff] = fminbnd(EqCondition,0.037,0.038);


% Corresponding stationary distribution
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

[mValueFunctionSS,mPolicyAssetSS,mPolicyConsSS,mIndexPolicyAssetSS,mPolicyLaborSS] = ...
    MultigridVFI_UBI(ttau,kkappa,llambda,rrho,rSS,ssigma,aalpha,A,depr,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,0,optAccelerator);

[mStationaryDistSS,expectAssetHoldingsSS] = StationaryDist(vGridAsset,nGridShock,mTransitionShock,mIndexPolicyAssetSS);
toc




%% Extended Path

T = 20;

vProd = ones(1,T);
vProd(1:10) = 10;

GuessR = rSS*ones(1,T);

bValueFunction1 = zeros(nGridAsset,nGridShock,T);
bValueFunction1(:,:,1) = mValueFunctionSS;
bValueFunction1(:,:,end) = mValueFunctionSS;

bPolicyAsset = zeros(nGridAsset,nGridShock,T);
bPolicyAsset(:,:,1) = mPolicyAssetSS;
bPolicyAsset(:,:,end) = mPolicyAssetSS;

bPolicyCons = zeros(nGridAsset,nGridShock,T);
bPolicyCons(:,:,1) = mPolicyConsSS;
bPolicyCons(:,:,end) = mPolicyConsSS;

bPolicyAssetIndex = zeros(nGridAsset,nGridShock,T);
bPolicyAssetIndex(:,:,1) = mIndexPolicyAssetSS;
bPolicyAssetIndex(:,:,end) = mIndexPolicyAssetSS;

bPolicyLabor = zeros(nGridAsset,nGridShock,T);
bPolicyLabor(:,:,1) = mPolicyLaborSS;
bPolicyLabor(:,:,end) = mPolicyLaborSS;


% Extended Path method

maxit = 1e04;
tol = 1e-6;
diff = 100;
it = 0;

while it<=maxit && diff>tol

    for t=1:T-2
    [bValueFunction1(:,:,T-t),bPolicyAsset(:,:,T-t),bPolicyCons(:,:,T-t),bPolicyAssetIndex(:,:,T-t),bPolicyLabor(:,:,T-t)] = ...
        VF_transition(bValueFunction1(:,:,T-t+1),ttau,llambda,kkappa,rrho,GuessR(T-t),aalpha,A,depr,ssigma,vGridAsset,vGridShock,mTransitionShock);
    end
end
    
    
% Backwards 
%for t=1:T-2
%    [bValueFunction1(:,:,T-t),bPolicyAsset(:,:,T-t),bPolicyCons(:,:,T-t),bPolicyAssetIndex(:,:,T-t),bPolicyLabor(:,:,T-t)] = ...
%        VF_transition(bValueFunction1(:,:,T-t+1),ttau,llambda,kkappa,rrho,GuessR(T-t),aalpha,A,depr,ssigma,vGridAsset,vGridShock,mTransitionShock);
%end




