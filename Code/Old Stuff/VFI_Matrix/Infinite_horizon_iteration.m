function [mValueFunction,mPolicyAsset,mPolicyCons,mIndexPolicyAsset] = ...
    Infinite_horizon_iteration(rrho,r,aalpha,A,depreciation,ssigma,vGridAsset,...
    vGridShock,mTransitionShock,mValueGuess)


% mValueGuess is the inital guess for the value function. Set to zero if
% there is no initial guess

wage = (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));

bbeta = 1/(1+rrho);
nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);

bAssetToday = repmat(reshape(vGridAsset,[nGridAsset,1]),[1,nGridShock,nGridAsset]);
bAssetNext = repmat(reshape(vGridAsset,[1,1,nGridAsset]),[nGridAsset,nGridShock,1]);
bShockToday = repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1,nGridAsset]);

bConsumption = wage*bShockToday + (1+r)*bAssetToday - bAssetNext;
%if sigma == 1
%    bConsumption(bConsumption<0) = 0;
%end
bConsumption(bConsumption<=0) = 0;

maxit = 1e04;
tol = 1e-5;
diff = 100;
it = 0;

if mValueGuess == 0
    mValueFunction = zeros(nGridAsset,nGridShock);
else
    mValueFunction = mValueGuess;
end

if ssigma == 1
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
    bUtility = ((bConsumption.^(1-ssigma))-1)/(1-ssigma);
    while it<=maxit && diff>tol
        bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
            [1,nGridShock,nGridAsset]),[nGridAsset,1,1]);
        bValue = (1-bbeta)*bUtility + bbeta*bContinuation;
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