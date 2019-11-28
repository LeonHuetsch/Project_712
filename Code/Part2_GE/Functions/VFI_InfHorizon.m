function [it,mValueFunction,mPolicyAsset,mPolicyCons] = ...
    VFI_InfHorizon(rrho,r,ssigma,aalpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,mValueGuess,optAccelerator)


wage = (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));
nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);
bbeta = 1/(1+rrho);

maxit = 1e04;
tol = 1e-06;
diff = 100;
it = 0;

if mValueGuess == 0
    mValueFunction = zeros(nGridAsset,nGridShock);
else
    mValueFunction = mValueGuess;
end
mPolicyAsset = zeros(nGridAsset,nGridShock);
mPolicyCons = zeros(nGridAsset,nGridShock);
mHelp = zeros(nGridAsset,nGridShock);

if ssigma == 1
    while it<=maxit && diff>tol
        it=it+1;
        mContValue = mValueFunction*mTransitionShock';

        for shockIndex=1:nGridShock
            shockToday = vGridShock(shockIndex);
            assetChoicePrev=1;

            for assetTodayIndex=1:nGridAsset
                assetToday = vGridAsset(assetTodayIndex);

                if mod(it,optAccelerator)==0 || it==1
                    valuePrev = -1000;
                    assetChoice = 1;

                    for assetNextIndex=assetChoicePrev:nGridAsset
                        assetNext = vGridAsset(assetNextIndex);
                        
                        cons = wage*shockToday + (1+r)*assetToday - assetNext;
                        if cons<=0.0 
                            value=-1000;
                        else
                            value=(1-bbeta)*log(cons) + bbeta*mContValue(assetNextIndex,shockIndex); 
                        end
                        if value>=valuePrev
                            valuePrev = value;
                            assetChoice = assetNextIndex;
                            assetChoicePrev = assetNextIndex;
                            consChoice = cons;
                        else      
                            break
                        end
                    end
                    mHelp(assetTodayIndex,shockIndex) = valuePrev;
                    mPolicyAsset(assetTodayIndex,shockIndex) = assetChoice; 
                    mPolicyCons(assetTodayIndex,shockIndex) = consChoice;
                else
                    mHelp(assetTodayIndex,shockIndex) = (1-bbeta)*log(mPolicyCons(assetTodayIndex,shockIndex))...
                        + bbeta*mContValue(mPolicyAsset(assetTodayIndex,shockIndex),shockIndex);
                end
            end
        end
        diff=max(max(abs(mHelp-mValueFunction)));
        mValueFunction = mHelp;        
    end		
    mPolicyAsset = vGridAsset(mPolicyAsset);
    
else 
    while it<=maxit && diff>tol
        it=it+1;
        mContValue = mValueFunction*mTransitionShock';

        for shockIndex=1:nGridShock
            shockToday = vGridShock(shockIndex);
            assetChoicePrev=1;

            for assetTodayIndex=1:nGridAsset
                assetToday = vGridAsset(assetTodayIndex);

                if mod(it,optAccelerator)==0 || it==1
                    valuePrev = -1000;
                    assetChoice = 1;
                    mHelp(assetTodayIndex,shockIndex)=-1e15;

                    for assetNextIndex=assetChoicePrev:nGridAsset
                        assetNext = vGridAsset(assetNextIndex);
                        cons = shockToday + (1+r)*assetToday - assetNext;
                        if cons<=0.0 
                            value=-1000;
                        else
                            value=(1-bbeta)*(cons^(1-ssigma)-1)/(1-ssigma) +...
                                bbeta*mContValue(assetNextIndex,shockIndex); 
                        end
                        if value>=valuePrev
                            valuePrev = value;
                            assetChoice = assetNextIndex;
                            assetChoicePrev = assetNextIndex;
                            consChoice = cons;
                        else      
                            break
                        end
                    end
                    mHelp(assetTodayIndex,shockIndex)=valuePrev;
                    mPolicyAsset(assetTodayIndex,shockIndex)=assetChoice; 
                    mPolicyCons(assetTodayIndex,shockIndex) = consChoice;
                else
                    mHelp(assetTodayIndex,shockIndex) = (1-bbeta)*log(mPolicyCons(assetTodayIndex,shockIndex))...
                        + bbeta*mContValue(mPolicyAsset(assetTodayIndex,shockIndex),shockIndex);
                end
            end
        end
        diff=max(max(abs(mHelp-mValueFunction)));
        mValueFunction = mHelp;        
    end	
    mPolicyAsset = vGridAsset(mPolicyAsset);
end
end