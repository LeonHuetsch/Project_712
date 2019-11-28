function [mValueFunction,mPolicyAsset,mPolicyCons] = ...
    Finite_horizon_iteration(rho,r,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt)


% Performs a value function iteration for finite horizon.
% If MortOpt=1, it includes conditional survival data from survs.txt 

if MortOpt == 0
    beta = 1/(1+rho);
    nGridAsset = length(vGridAsset);
    nGridShock = length(vGridShock);

    bAssetToday = repmat(reshape(vGridAsset,[nGridAsset,1]),[1,nGridShock,nGridAsset]);
    bAssetNext = repmat(reshape(vGridAsset,[1,1,nGridAsset]),[nGridAsset,nGridShock,1]);
    bShockToday = repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1,nGridAsset]);

    bConsumption = bShockToday + (1+r)*bAssetToday - bAssetNext;
    bConsumption(bConsumption<0) = 0;

    mValueFunction = zeros(nGridAsset,nGridShock,nPeriod);
    mPolicyAsset = zeros(nGridAsset,nGridShock,nPeriod);
    mPolicyCons = zeros(nGridAsset,nGridShock,nPeriod);


    for t=2:nPeriod
        bContinuation = repmat(reshape(mTransitionShock*mValueFunction(:,:,t-1)',...
            [1,nGridShock,nGridAsset]),[nGridAsset,1,1]);
        bValue = (1-beta)*log(bConsumption) + beta*bContinuation;
        bValue(bConsumption<0) = -1e20;

        [mValueFunction(:,:,t),indexMax] = max(bValue,[],3); 
        mPolicyAsset(:,:,t) = vGridAsset(indexMax);
        mPolicyCons(:,:,t) = bShockToday(:,:,1) + (1+r)*bAssetToday(:,:,1)...
            - mPolicyAsset(:,:,t);
        %mPolicyCons(:,:,T) = bConsumption(:,:,indexMax);
    end

    
else
    
    inpath = 'Data/';
    vMortData = load([inpath,'survs.txt']); 
    
    vSurvivalProb = vMortData(end:-1:1);    
    beta = 1/(1+rho);
    vEffectiveDiscout = beta*vSurvivalProb;
    
    nGridAsset = length(vGridAsset);
    nGridShock = length(vGridShock);

    bAssetToday = repmat(reshape(vGridAsset,[nGridAsset,1]),[1,nGridShock,nGridAsset]);
    bAssetNext = repmat(reshape(vGridAsset,[1,1,nGridAsset]),[nGridAsset,nGridShock,1]);
    bShockToday = repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1,nGridAsset]);

    bConsumption = bShockToday + (1+r)*bAssetToday - bAssetNext;
    bConsumption(bConsumption<0) = 0;

    mValueFunction = zeros(nGridAsset,nGridShock,nPeriod);
    mPolicyAsset = zeros(nGridAsset,nGridShock,nPeriod);
    mPolicyCons = zeros(nGridAsset,nGridShock,nPeriod);


    for t=2:nPeriod
        bContinuation = repmat(reshape(mTransitionShock*mValueFunction(:,:,t-1)',...
            [1,nGridShock,nGridAsset]),[nGridAsset,1,1]);
        bValue = log(bConsumption) + ...
            vEffectiveDiscout(t-1)*bContinuation;
        bValue(bConsumption<0) = -1e20;

        [mValueFunction(:,:,t),indexMax] = max(bValue,[],3); 
        mPolicyAsset(:,:,t) = vGridAsset(indexMax);
        mPolicyCons(:,:,t) = bShockToday(:,:,1) + (1+r)*bAssetToday(:,:,1)...
            - mPolicyAsset(:,:,t);
        %mPolicyCons(:,:,T) = bConsumption(:,:,indexMax);
    end
    
end

end