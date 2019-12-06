function [mValueFunctionToday,mPolicyAsset,mPolicyCons,mPolicyAssetIndex] = ...
    VF_transition(mValueFunctionNext,rrho,r,aalpha,A,depr,ssigma,vGridAsset,vGridShock,mTransitionShock)


wage = (1-aalpha)*A*(aalpha*A/(r+depr))^(aalpha/(1-aalpha));

bbeta = 1/(1+rrho);
nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);

bAssetToday = repmat(reshape(vGridAsset,[nGridAsset,1]),...
    [1,nGridShock,nGridAsset]);
bShockToday = repmat(reshape(vGridShock,[1,nGridShock]),...
    [nGridAsset,1,nGridAsset]);
bAssetNext = repmat(reshape(vGridAsset,[1,1,nGridAsset]),...
    [nGridAsset,nGridShock]);
bConsumption = wage*bShockToday + (1+r)*bAssetToday...
    - bAssetNext;

bConsumption(bConsumption<=0) = 0;

if ssigma == 1
    bUtility = log(bConsumption);
    bContinuation = repmat(reshape(mTransitionShock*mValueFunctionNext',...
        [1,nGridShock,nGridAsset]),[nGridAsset,1,1]);
    bValue = bUtility + bbeta*bContinuation;
    bValue(bConsumption<=0) = -1e20;

    mValueFunctionToday=max(bValue,[],3);
else
    bUtility = ((bConsumption.^(1-ssigma))-1)/(1-ssigma);
    bContinuation = repmat(reshape(mTransitionShock*mValueFunctionNext',...
        [1,nGridShock,nGridAsset]),[nGridAsset,1,1]);
    bValue = bUtility + bbeta*bContinuation;
    bValue(bConsumption<0) = -1e20;

    mValueFunctionToday=max(bValue,[],3);
end

[~,mPolicyIndex]=max(bValue,[],3);

mPolicyAssetIndex = rem(mPolicyIndex,nGridAsset);
mPolicyAssetIndex(mPolicyAssetIndex==0)=nGridAsset;
mPolicyAsset=vGridAsset(mPolicyAssetIndex);

mPolicyCons = wage*bShockToday(:,:,1) + (1+r)*bAssetToday(:,:,1)...
        - mPolicyAsset;

end