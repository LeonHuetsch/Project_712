function [mValueFunction,mPolicyAsset,mPolicyCons] = ...
    Finite_horizon_iteration(r,vGridAsset,vGridShock,mTransitionShock,nPeriod)


nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);

mValueFunction = zeros(nGridAsset,nGridShock,nPeriod);
mPolicyAsset = zeros(nGridAsset,nGridShock,nPeriod);
mPolicyCons = zeros(nGridAsset,nGridShock,nPeriod);
mHelp = zeros(nPeriod,nGridAsset,nGridShock);

for T=1:nPeriod
        for shockToday=1:nGridShock
            capitalChoice=1;
        for assetToday=1:nGridAsset
            vHelpPrev=-1e12;
            mHelp(assetToday,shockToday,T)=-1e15;
            for assetNext=capitalChoice:nGridAsset
                cons = vGridShock(shockToday) + (1+r)*vGridAsset(assetToday)...
                    - vGridAsset(assetNext);
                if cons<=0.0 
                    vhelp=-10^12;
                else
                    vhelp=log(cons) +...
                        beta*sum(mValueFunction(assetNext,:).*mTransitionShock(shockToday,:)); 
                end
                if vhelp<vHelpPrev
                    break
                end
                if vhelp>=mHelp(assetToday,shockToday,T) 
                    mHelp(assetToday,shockToday,T)=vhelp;
                    mPolicyAsset(assetToday,shockToday,T)=vGridAsset(assetNext); 
                    mPolicyCons(assetToday,shockToday,T)=cons;
                    capitalChoice=assetNext;
                end
                vHelpPrev=vhelp;
            end
        end
        end
        mValueFunction(assetToday,shockToday,T) = mHelp(assetToday,shockToday,T);        
end		 

end