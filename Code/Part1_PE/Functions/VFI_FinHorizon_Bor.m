function [mValueFunction,mPolicyAsset,mPolicyCons] = ...
    VFI_FinHorizon(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,nPeriod,mortOpt)


nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);

minAsset = vGridAsset(1);
maxAsset = vGridAsset(end);

bbeta = 1/(1+rrho);

mValueFunction = zeros(nGridAsset,nGridShock,nPeriod);
mPolicyAsset = zeros(nGridAsset,nGridShock,nPeriod);
mPolicyCons = zeros(nGridAsset,nGridShock,nPeriod);


if mortOpt == 0
    if ssigma==1
        for t=2:nPeriod
            mContValue = mValueFunction(:,:,t-1)*mTransitionShock';

            for shockIndex=1:nGridShock
                shockToday = vGridShock(shockIndex);
                assetChoice = minAsset;

                vContValue = griddedInterpolant(vGridAsset',mContValue(:,shockIndex));
                for assetTodayIndex=1:nGridAsset
                    assetToday = vGridAsset(assetTodayIndex);
                    
                    if t==2
                        cons = shockToday+(1+r)*assetToday;
                        assetChoice = 0;
                        if cons<=0
                            value = -10000;
                        else 
                            value = log(cons);
                        end
                        mValueFunction(assetTodayIndex,shockIndex,t) = value;
                        mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                        mPolicyCons(assetTodayIndex,shockIndex,t) = cons;   
                        
                    else
                    
                        maxi = shockToday + (1+r)*assetToday;       % ensures positive consumption

                        valuefun = @(at) -1*(log(shockToday+(1+r)*assetToday-at) + bbeta*vContValue(at));
                        [assetChoice,value] = fminbnd(valuefun,assetChoice,maxi);%,options);

                        mValueFunction(assetTodayIndex,shockIndex,t) = -value;
                        mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                        mPolicyCons(assetTodayIndex,shockIndex,t) = shockToday + (1+r)*assetToday - assetChoice;
                    end
                end
            end        
        end	
    else
        for t=2:nPeriod
            mContValue = mValueFunction(:,:,t-1)*mTransitionShock';

            for shockIndex=1:nGridShock
                shockToday = vGridShock(shockIndex);
                assetChoice = minAsset;

                vContValue = griddedInterpolant(vGridAsset',mContValue(:,shockIndex));
                for assetTodayIndex=1:nGridAsset
                    assetToday = vGridAsset(assetTodayIndex);
                    
                    if t==2
                        cons = shockToday+(1+r)*assetToday;
                        assetChoice = 0;
                        if cons<=0
                            value = -1000;
                        else 
                            value = (cons^(1-ssigma)-1)/(1-ssigma);
                        end
                        mValueFunction(assetTodayIndex,shockIndex,t) = value;
                        mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                        mPolicyCons(assetTodayIndex,shockIndex,t) = cons;   
                        
                    else
                        maxi = shockToday + (1+r)*assetToday;

                        valuefun = @(at) -1*((((shockToday+(1+r)*assetToday-at)^(1-ssigma)-1)/(1-ssigma)) + bbeta*vContValue(at));

                        [assetChoice,value] = fminbnd(valuefun,assetChoice,maxi);%,options);

                        mValueFunction(assetTodayIndex,shockIndex,t) = -value;
                        mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                        mPolicyCons(assetTodayIndex,shockIndex,t) = shockToday + (1+r)*assetToday - assetChoice;
                    end
                end
            end        
        end	
    end
  
else
    if ssigma==1
        inpath = 'Data/';
        vMortData = load([inpath,'survs.txt']); 

        vSurvivalProb = vMortData(end:-1:1);    
        vEffectiveDiscout = bbeta*vSurvivalProb;

        for t=2:nPeriod
            mContValue = mValueFunction(:,:,t-1)*mTransitionShock';

            for shockIndex=1:nGridShock
                shockToday = vGridShock(shockIndex);
                assetChoice = minAsset;
                
                vContValue = griddedInterpolant(vGridAsset',mContValue(:,shockIndex));
                for assetTodayIndex=1:nGridAsset
                    assetToday = vGridAsset(assetTodayIndex);
                    
                    if t==2
                        cons = shockToday+(1+r)*assetToday;
                        assetChoice = 0;
                        if cons<=0
                            value = -1000;
                        else 
                            value = log(cons);
                        end
                        mValueFunction(assetTodayIndex,shockIndex,t) = value;
                        mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                        mPolicyCons(assetTodayIndex,shockIndex,t) = cons;   
                        
                    else
                        maxi = shockToday + (1+r)*assetToday;

                        valuefun = @(at) -1*(log(shockToday+(1+r)*assetToday-at) + vEffectiveDiscout(t-1)*vContValue(at));

                        [assetChoice,value] = fminbnd(valuefun,assetChoice,maxi);%,options);

                        mValueFunction(assetTodayIndex,shockIndex,t) = -value;
                        mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                        mPolicyCons(assetTodayIndex,shockIndex,t) = shockToday + (1+r)*assetToday - assetChoice;
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
                assetChoice = minAsset;
                
                vContValue = griddedInterpolant(vGridAsset',mContValue(:,shockIndex));
                for assetTodayIndex=1:nGridAsset
                    assetToday = vGridAsset(assetTodayIndex);
                    
                    if t==2
                        cons = shockToday+(1+r)*assetToday;
                        assetChoice = 0;
                        if cons<=0
                            value = -1000;
                        else 
                            value = (cons^(1-ssigma)-1)/(1-ssigma);
                        end
                        mValueFunction(assetTodayIndex,shockIndex,t) = value;
                        mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                        mPolicyCons(assetTodayIndex,shockIndex,t) = cons;   
                        
                    else
                    maxi = shockToday + (1+r)*assetToday;
                    
                    valuefun = @(at) -1*((((shockToday+(1+r)*assetToday-at)^(1-ssigma)-1)/(1-ssigma)) + vEffectiveDiscout(t-1)*vContValue(at));
                    
                    [assetChoice,value] = fminbnd(valuefun,assetChoice,maxi);%,options);

                    mValueFunction(assetTodayIndex,shockIndex,t) = -value;
                    mPolicyAsset(assetTodayIndex,shockIndex,t) = assetChoice; 
                    mPolicyCons(assetTodayIndex,shockIndex,t) = shockToday + (1+r)*assetToday - assetChoice;
                    end
                end
            end        
        end	
    end
end
end