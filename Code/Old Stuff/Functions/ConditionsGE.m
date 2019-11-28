function vEqConditions = ConditionsGE(nGridAsset,minGridAsset,...
    maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt,...
    rho,alpha,A,depreciation,sigma,mValueGuess,lambda,r,kappa,tau) 

%if length(EqParameters)==1
%    r = EqParameters;
%    kappa = 0;
%    tau = 0;    
%else
    %r = EqParameters(1);
    %kappa = EqParameters(2);
    %tau = EqParameters(3);
    %tau=0;
    %tau=0;
%end

% Given interest rate r, tax rate tau and disutility of working kappa
% the function calculates savings, tax revenue and the share of population
% working.
% It returns a vector with the difference between these variables and capital,
% demand, lambda and the required working share 0.8, respectively. 
% Fsolve can minimize this vector to find the equilibrium r, lambda, kappa

wage = (1-alpha)*A*(alpha*A/(r+depreciation))^(alpha/(1-alpha));

%% Savings

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,...
    logShockAverage,truncOpt);


[~,mPolicyAsset,~,~,mPolicyLabor] = ...
    VFiteration_UBI(tau,lambda,kappa,rho,r,alpha,A,depreciation,sigma,vGridAsset,...
    vGridShock,mTransitionShock,mValueGuess);


[mStationaryDist,expectAssetHoldings] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset);


%% Capital demand 

effectiveLaborSup = sum(sum(mStationaryDist.*mPolicyLabor.*...
    repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])));
capitalDemand = (alpha*A/(r+depreciation))^(1/(1-alpha))*effectiveLaborSup;

%% Share of population working

workingShare = sum(sum(mStationaryDist.*mPolicyLabor));

%% Tax revenue

govtRevenue = tau*wage*(sum(sum(mStationaryDist.*(mPolicyLabor.*...
    repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])))));

%% Return vector of equilibrium conditions

%vEqConditions = zeros(1,2);
%if length(EqParameters)==1
%    vEqConditions = capitalDemand - expectAssetHoldings;
%else
    %vEqConditions(2) = capitalDemand - expectAssetHoldings;
    %vEqConditions = (workingShare - 0.8)^2 + (capitalDemand - expectAssetHoldings)^2;
vEqConditions(1) = capitalDemand - expectAssetHoldings;
vEqConditions(2) = 35*(workingShare - 0.8);
vEqConditions(3) = 30*(govtRevenue - lambda);
%end

end