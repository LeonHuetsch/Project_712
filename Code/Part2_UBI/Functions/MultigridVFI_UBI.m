function [mValueFunction,mPolicyAsset,mPolicyCons,mIndexPolicyAsset,mPolicyLabor]...
    = MultigridVFI_UBI(ttau,kkappa,llambda,rrho,r,ssigma,aalpha,A,depreciation,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,mGuessVF,optAccelerator)


nGridShock = length(vGridShock);

% Interpolate initial guess so it fits first gridsize for assets
if mGuessVF ~= 0
    F = griddedInterpolant(mGuessVF);
    mRows = repmat(reshape(linspace(1,size(mGuessVF,1),vMultiSteps(1)),[vMultiSteps(1),1]),[1,nGridShock]);
    mColumns = repmat(1:nGridShock,[vMultiSteps(1),1]);
    mGuessVF = F(mRows,mColumns);
end
for step=1:length(vMultiSteps)
    if step==1
        vGridAsset = linspace(minAsset,maxAsset,vMultiSteps(step));
        nGridAsset = length(vGridAsset);
        [mValueFunction,mPolicyAsset,mPolicyCons,mIndexPolicyAsset,mPolicyLabor] = VFI_UBI(ttau,kkappa,llambda,rrho,r,ssigma,aalpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,mGuessVF,optAccelerator);
    else
        F = griddedInterpolant(mValueFunction);
        mFinerRows = repmat(reshape(linspace(1,nGridAsset,vMultiSteps(step)),[vMultiSteps(step),1]),[1,nGridShock]);
        mColumns = repmat(1:nGridShock,[vMultiSteps(step),1]);
        mInitialGuess = F(mFinerRows,mColumns);
        
        vGridAsset = linspace(minAsset,maxAsset,vMultiSteps(step));
        nGridAsset = length(vGridAsset);

        [mValueFunction,mPolicyAsset,mPolicyCons,mIndexPolicyAsset,mPolicyLabor] = VFI_UBI(ttau,kkappa,llambda,rrho,r,ssigma,aalpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,mInitialGuess,optAccelerator);
    end
end
end