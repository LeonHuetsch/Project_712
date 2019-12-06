function [it,mValueFunction,mPolicyAsset,mPolicyCons]...
    = MultigridVFI_InfHorizon(rrho,r,ssigma,aalpha,A,depreciation,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,mGuessVF,optAccelerator)

nGridShock = length(vGridShock);
for step=1:length(vMultiSteps)
    if step==1
        vGridAsset = linspace(minAsset,maxAsset,vMultiSteps(step));
        nGridAsset = length(vGridAsset);
        [it,mValueFunction,mPolicyAsset,mPolicyCons] = VFI_InfHorizon(rrho,r,ssigma,aalpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,mGuessVF,optAccelerator);
    else
        F = griddedInterpolant(mValueFunction);
        mFinerRows = repmat(reshape(linspace(1,nGridAsset,vMultiSteps(step)),[vMultiSteps(step),1]),[1,nGridShock]);
        mFinerColumns = repmat(1:nGridShock,[vMultiSteps(step),1]);
        mInitialGuess = F(mFinerRows,mFinerColumns);
        
        vGridAsset = linspace(minAsset,maxAsset,vMultiSteps(step));
        nGridAsset = length(vGridAsset);

        [~,mValueFunction,mPolicyAsset,mPolicyCons] = VFI_InfHorizon(rrho,r,ssigma,aalpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,mInitialGuess,optAccelerator);
    end
end
end