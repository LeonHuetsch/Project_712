function [mValueFunction,mPolicyAsset,mPolicyCons,mPolicyAssetIndex,mPolicyLabor] = VFI_UBI(ttau,kkappa,llambda,rrho,r,ssigma,alpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,mValueGuess,optAccelerator)


% mValueGuess is the inital guess for the value function. Set to zero if
% there is no initial guess. 
nGridLabor = 2;
wage = (1-alpha)*A*(alpha*A/(r+depreciation))^(alpha/(1-alpha));

bbeta = 1/(1+rrho);
nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);

if mValueGuess == 0
    mValueFunction = zeros(nGridAsset,nGridShock);
else
    mValueFunction = mValueGuess;
end

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


maxit = 1e04;
tol = 1e-7;
diff = 100;
it = 0;

if ssigma == 1
    bUtility = log(bConsumption) - kkappa*bLaborToday;
    while it<=maxit && diff>tol
        it=it+1;
        if mod(it,optAccelerator)==0 || it==1
            bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
                [1,nGridShock,nGridAsset]),[nGridAsset,1,nGridLabor]);
            bValue = bUtility + bbeta*bContinuation;
            bValue(bConsumption<=0) = -1e20;

            [mHelp,mPolicyIndex] = max(bValue,[],3);

            %diff=max(max(abs(mHelp-mValueFunction)));    
            %mValueFunction=mHelp;   
            
            mPolicyAssetIndex = rem(mPolicyIndex,nGridAsset);
            mPolicyAssetIndex(mPolicyAssetIndex==0)=nGridAsset;
            mPolicyAsset=vGridAsset(mPolicyAssetIndex);
            mPolicyLabor=mPolicyIndex<=nGridAsset;
            mPolicyCons = (1-ttau)*wage*bShockToday(:,:,1).*mPolicyLabor + (1+r)*bAssetToday(:,:,1)...
                - mPolicyAsset + llambda;
        else
            mContValue = mValueFunction*mTransitionShock'; 
            for shockIndex=1:nGridShock
                for assetTodayIndex=1:nGridAsset
                    
                    mHelp(assetTodayIndex,shockIndex) = (log(mPolicyCons(assetTodayIndex,shockIndex))-kkappa*mPolicyLabor(assetTodayIndex,shockIndex))...
                        + bbeta*mContValue(mPolicyAssetIndex(assetTodayIndex,shockIndex),shockIndex);          
                end
            end
        end
        diff=max(max(abs(mHelp-mValueFunction)));    
        mValueFunction=mHelp;  
    end	

else
    bUtility = ((bConsumption.^(1-ssigma))-1)/(1-ssigma) - kkappa*bLaborToday;
    while it<=maxit && diff>tol
        bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
            [1,nGridShock,nGridAsset]),[nGridAsset,1,nGridLabor]);
        bValue = bUtility + bbeta*bContinuation;
        bValue(bConsumption<0) = -1e20;

        mHelp=max(bValue,[],3);

        diff=max(max(abs(mHelp-mValueFunction)));    
        mValueFunction=mHelp;   
        it=it+1;
    end	    
end

[~,mPolicyIndex]=max(bValue,[],3);

mPolicyAssetIndex = rem(mPolicyIndex,nGridAsset);
mPolicyAssetIndex(mPolicyAssetIndex==0)=nGridAsset;
mPolicyAsset=vGridAsset(mPolicyAssetIndex);

mPolicyLabor=mPolicyIndex<=nGridAsset;

mPolicyCons = (1-ttau)*wage*bShockToday(:,:,1).*mPolicyLabor + (1+r)*bAssetToday(:,:,1)...
        - mPolicyAsset + llambda;

end