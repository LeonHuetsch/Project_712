function [mValueFunctionToday,mPolicyAsset,mPolicyCons,mPolicyAssetIndex,mPolicyLabor] = ...
    VF_transition(mValueFunctionNext,ttau,llambda,kkappa,rrho,r,aalpha,A,depr,ssigma,vGridAsset,vGridShock,mTransitionShock)


nGridLabor = 2;
wage = (1-aalpha)*A*(aalpha*A/(r+depr))^(aalpha/(1-aalpha));

bbeta = 1/(1+rrho);
nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);

bAssetToday = repmat(reshape(vGridAsset,[nGridAsset,1]),...
    [1,nGridShock,nGridAsset*nGridLabor]);
bShockToday = repmat(reshape(vGridShock,[1,nGridShock]),...
    [nGridAsset,1,nGridAsset*nGridLabor]);
bAssetNext = repmat(reshape(vGridAsset,[1,1,nGridAsset]),...
    [nGridAsset,nGridShock,nGridLabor]);
bLaborToday = zeros(nGridAsset,nGridShock,nGridAsset*nGridLabor);
bLaborToday(:,:,1:nGridAsset) = 1;    

bConsumption = (1-ttau)*wage*bShockToday.*bLaborToday + (1+r)*bAssetToday...
    - bAssetNext + llambda;

bConsumption(bConsumption<=0) = 0;

if ssigma == 1
    bUtility = (1-bbeta)*(log(bConsumption) - kkappa*bLaborToday);
    bContinuation = repmat(reshape(mTransitionShock*mValueFunctionNext',...
        [1,nGridShock,nGridAsset]),[nGridAsset,1,nGridLabor]);
    bValue = bUtility + bbeta*bContinuation;
    bValue(bConsumption<=0) = -1e20;

    mValueFunctionToday=max(bValue,[],3);
else
    bUtility = (1-bbeta)*(((bConsumption.^(1-ssigma))-1)/(1-ssigma) - kkappa*bLaborToday);
    bContinuation = repmat(reshape(mTransitionShock*mValueFunctionNext',...
        [1,nGridShock,nGridAsset]),[nGridAsset,1,nGridLabor]);
    bValue = bUtility + bbeta*bContinuation;
    bValue(bConsumption<0) = -1e20;

    mValueFunctionToday=max(bValue,[],3);
end

[~,mPolicyIndex]=max(bValue,[],3);

mPolicyAssetIndex = rem(mPolicyIndex,nGridAsset);
mPolicyAssetIndex(mPolicyAssetIndex==0)=nGridAsset;
mPolicyAsset=vGridAsset(mPolicyAssetIndex);

mPolicyLabor=mPolicyIndex<=nGridAsset;

mPolicyCons = (1-ttau)*wage*bShockToday(:,:,1).*mPolicyLabor + (1+r)*bAssetToday(:,:,1)...
        - mPolicyAsset + llambda;

end