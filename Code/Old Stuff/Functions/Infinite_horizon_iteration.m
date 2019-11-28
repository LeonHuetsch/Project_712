function [mValueFunction,mPolicyAsset,mPolicyCons,mIndexPolicyAsset] = ...
    Infinite_horizon_iteration(rho,r,aalpha,A,depreciation,sigma,vGridAsset,...
    vGridShock,mTransitionShock,mValueGuess)


% mValueGuess is the inital guess for the value function. Set to zero if
% there is no initial guess

wage = (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));
wage = 1;

bbeta = 1/(1+rho);
nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);

bAssetToday = repmat(reshape(vGridAsset,[nGridAsset,1]),[1,nGridShock,nGridAsset]);
bAssetNext = repmat(reshape(vGridAsset,[1,1,nGridAsset]),[nGridAsset,nGridShock,1]);
bShockToday = repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1,nGridAsset]);

bConsumption = wage*bShockToday + (1+r)*bAssetToday - bAssetNext;

bConsumption(bConsumption<=0) = 0;

maxit = 1e04;
tol = 1e-06;
diff = 100;
it = 0;

if mValueGuess == 0
    mValueFunction = zeros(nGridAsset,nGridShock);
else
    mValueFunction = mValueGuess;
end

if sigma == 1
    bUtility = log(bConsumption);
    while it<=maxit && diff>tol
        bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
            [1,nGridShock,nGridAsset]),[nGridAsset,1,1]);
        bValue = (1-bbeta)*bUtility + bbeta*bContinuation;
        bValue(bConsumption<=0) = -1e20;

        mHelp=max(bValue,[],3);

        diff=max(max(abs(mHelp-mValueFunction)));    
        mValueFunction=mHelp;   
        it=it+1;
    end	
    
else
    bUtility = ((bConsumption.^(1-sigma))-1)/(1-sigma);
    while it<=maxit && diff>tol
        bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
            [1,nGridShock,nGridAsset]),[nGridAsset,1,1]);
        bValue = bUtility + bbeta*bContinuation;
        bValue(bConsumption<0) = -1e20;

        mHelp=max(bValue,[],3);

        diff=max(max(abs(mHelp-mValueFunction)));    
        mValueFunction=mHelp;   
        it=it+1;
    end	    
end

[~,mIndexPolicyAsset]=max(bValue,[],3);
mPolicyAsset=vGridAsset(mIndexPolicyAsset);
mPolicyCons = wage*bShockToday(:,:,1) + (1+r)*bAssetToday(:,:,1)...
        - mPolicyAsset;

end