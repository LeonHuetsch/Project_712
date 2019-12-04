function EqCondition = ConditionsGE(nGridAsset,minAsset,maxAsset,vMultiSteps,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt,rrho,aalpha,A,depreciation,ssigma,mValueGuess,llambda,r,kkappa,ttau,optAccelerator)


% Given interest rate r, tax rate tau and disutility of working kappa
% the function calculates savings, tax revenue and the share of population
% working.
% It returns a vector with the difference between these variables and capital,
% demand, lambda and the required working share 0.8, respectively. 
% Fsolve can minimize this vector to find the equilibrium r, lambda, kappa

%wage = (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));

%% Savings
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

[~,~,~,mIndexPolicyAsset,mPolicyLabor] = MultigridVFI_UBI(ttau,kkappa,llambda,rrho,r,ssigma,aalpha,A,depreciation,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,mValueGuess,optAccelerator);

[mStationaryDist,expectAssetHoldings] = StationaryDist(vGridAsset,nGridShock,mTransitionShock,mIndexPolicyAsset);


%% Capital demand 
effectiveLaborSup = sum(sum(mStationaryDist.*mPolicyLabor.*repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])));
capitalDemand = (aalpha*A/(r+depreciation))^(1/(1-aalpha))*effectiveLaborSup;


%% Share of population working
%workingShare = sum(sum(mStationaryDist.*mPolicyLabor));


%% Tax revenue
%govtRevenue = ttau*wage*(sum(sum(mStationaryDist.*(mPolicyLabor.*...
%    repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])))));


%% Return vector of equilibrium conditions
%{
if llambda == 0
    vEqConditions(1) = capitalDemand - expectAssetHoldings;
    vEqConditions(2) = 20*(workingShare - 0.8);
else
    vEqConditions(1) = capitalDemand - expectAssetHoldings;
    vEqConditions(2) = 20*(govtRevenue - llambda);
end
%}

%% Return vector of equilibrium conditions
EqCondition = capitalDemand - expectAssetHoldings;



end