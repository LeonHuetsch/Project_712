function [mValueFunction,mPolicyAsset,mPolicyCons] = ...
    VFI_FinHorizon(rrho,r,vGridAsset,vGridShock,mTransitionShock,nPeriod,mortOpt)


nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);
bbeta = 1/(1+rrho);

mValueFunction = zeros(nGridAsset,nGridShock,nPeriod);
mPolicyAsset = zeros(nGridAsset,nGridShock,nPeriod);
mPolicyCons = zeros(nGridAsset,nGridShock,nPeriod);


if mortOpt == 0
    for t=2:nPeriod
        mContValue = mValueFunction(:,:,t-1)*mTransitionShock';
        
        for shockIndex=1:nGridShock
            shockToday = vGridShock(shockIndex);
            assetChoicePrev=1;

            for assetTodayIndex=1:nGridAsset
                assetToday = vGridAsset(assetTodayIndex);
                valuePrev = -1000;
                assetChoice = 1;

                for assetNextIndex=assetChoicePrev:nGridAsset
                    assetNext = vGridAsset(assetNextIndex);

                    cons = shockToday + (1+r)*assetToday - assetNext;
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
                    mValueFunction(assetTodayIndex,shockIndex,t) = valuePrev;
                    mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                    mPolicyCons(assetTodayIndex,shockIndex,t) = consChoice;
                end
            end
        end        
    end	
    
else
    inpath = 'Data/';
    vMortData = load([inpath,'survs.txt']); 
    
    vSurvivalProb = vMortData(end:-1:1);    
    vEffectiveDiscout = bbeta*vSurvivalProb;
    
    for t=2:nPeriod
        mContValue = mValueFunction(:,:,t-1)*mTransitionShock';
        
        for shockIndex=1:nGridShock
            shockToday = vGridShock(shockIndex);
            assetChoicePrev=1;

            for assetTodayIndex=1:nGridAsset
                assetToday = vGridAsset(assetTodayIndex);
                valuePrev = -1000;
                assetChoice = 1;

                for assetNextIndex=assetChoicePrev:nGridAsset
                    assetNext = vGridAsset(assetNextIndex);

                    cons = shockToday + (1+r)*assetToday - assetNext;
                    if cons<=0.0 
                        value=-1000;
                    else
                        value=(1-bbeta)*log(cons) + vEffectiveDiscout(t-1)*mContValue(assetNextIndex,shockIndex); 
                    end
                    if value>=valuePrev
                        valuePrev = value;
                        assetChoice = assetNextIndex;
                        assetChoicePrev = assetNextIndex;
                        consChoice = cons;
                    else      
                        break
                    end
                    mValueFunction(assetTodayIndex,shockIndex,t) = valuePrev;
                    mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                    mPolicyCons(assetTodayIndex,shockIndex,t) = consChoice;
                end
            end
        end        
    end	
end