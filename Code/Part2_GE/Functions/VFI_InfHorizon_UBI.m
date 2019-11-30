function [it,mValueFunction1,mPolicyAsset1,mPolicyCons1,mPolicyLabor1] = ...
    VFI_InfHorizon_UBI(tau,kappa,lambda,rrho,r,ssigma,aalpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,mValueGuess,optAccelerator)


wage = (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));
nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);
bbeta = 1/(1+rrho);

maxit = 1e4;
tol = 1e-06;
diff = 100;
it = 0;

if mValueGuess == 0
    mValueFunction1 = zeros(nGridAsset,nGridShock);
else
    mValueFunction1 = mValueGuess;
end
mPolicyAsset1 = zeros(nGridAsset,nGridShock);
mPolicyCons1 = zeros(nGridAsset,nGridShock);
mPolicyLabor1 = zeros(nGridAsset,nGridShock);
mHelp = zeros(nGridAsset,nGridShock);

if ssigma == 1
    while it<=maxit && diff>tol
        it=it+1;
        mContValue = mValueFunction1*mTransitionShock';

        for shockIndex=1:nGridShock
            shockToday = vGridShock(shockIndex);
            assetChoicePrev=1;
            upper = nGridAsset;

            for assetTodayIndex=1:nGridAsset
                assetToday = vGridAsset(assetTodayIndex);

                if mod(it,optAccelerator)==0 || it==1
                    valuePrev = -1e20;
                    assetChoice = 1;
                    consChoice = 0;

                    for assetNextIndex=assetChoicePrev:upper
                        assetNext = vGridAsset(assetNextIndex);
                        
                        consWork = (1-tau)*wage*shockToday+(1+r)*assetToday+lambda-assetNext;
                        consLeisure =  (1+r)*assetToday+lambda-assetNext;
                        if consWork<=0
                            break
                        end
                        
                        consWork = max(0,consWork);
                        consLeisure =  max(0,consLeisure);
                        
                        workVal = log(consWork) - kappa;
                        leisureVal = log(consLeisure);
                        
                        if workVal>=leisureVal
                            labor = 1;
                            cons = consWork;
                            value = workVal + bbeta*mContValue(assetNextIndex,shockIndex); 
                        else
                            labor = 0;
                            cons = consLeisure; 
                            value = leisureVal + bbeta*mContValue(assetNextIndex,shockIndex); 
                        end

                        if value>valuePrev
                            valuePrev = value;
                            assetChoice = assetNextIndex;
                            assetChoicePrev = max(assetNextIndex-10,1);
                            upper = min(assetNextIndex+10,nGridAsset);
                            consChoice = cons;
                            laborChoice = labor;
                        end
                    end
                    mHelp(assetTodayIndex,shockIndex) = valuePrev;
                    mPolicyAsset1(assetTodayIndex,shockIndex) = assetChoice; 
                    mPolicyCons1(assetTodayIndex,shockIndex) = consChoice;
                    mPolicyLabor1(assetTodayIndex,shockIndex) = laborChoice;
                else
                    mHelp(assetTodayIndex,shockIndex) = (log(mPolicyCons1(assetTodayIndex,shockIndex))-kappa*mPolicyLabor1(assetTodayIndex,shockIndex))...
                        + bbeta*mContValue(mPolicyAsset1(assetTodayIndex,shockIndex),shockIndex);            
                end
            end
        end
        diff=max(max(abs(mHelp-mValueFunction1)));
        mValueFunction1 = mHelp;        
    end		
    %mPolicyAsset = vGridAsset(mPolicyAsset);
    
else 
    while it<=maxit && diff>tol
        it=it+1;
        mContValue = mValueFunction1*mTransitionShock';

        for shockIndex=1:nGridShock
            shockToday = vGridShock(shockIndex);
            assetChoicePrev=1;

            for assetTodayIndex=1:nGridAsset
                assetToday = vGridAsset(assetTodayIndex);

                if mod(it,optAccelerator)==0 || it==1
                    valuePrev = -1e06;
                    assetChoice = 1;
                    consChoice = 0;

                    for assetNextIndex=assetChoicePrev:nGridAsset
                        assetNext = vGridAsset(assetNextIndex);
                        
                        cons = wage*shockToday + (1+r)*assetToday - assetNext;
                        if cons<=0.0 
                            value=-1e07;
                        else
                            value=(1-bbeta)*(cons^(1-ssigma)-1)/(1-ssigma) + bbeta*mContValue(assetNextIndex,shockIndex); 
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
                    mPolicyAsset1(assetTodayIndex,shockIndex)=assetChoice; 
                    mPolicyCons1(assetTodayIndex,shockIndex) = consChoice;
                else
                    mHelp(assetTodayIndex,shockIndex) = (1-bbeta)*(mPolicyCons1(assetTodayIndex,shockIndex)^(1-ssigma)-1)/(1-ssigma)...
                        + bbeta*mContValue(mPolicyAsset1(assetTodayIndex,shockIndex),shockIndex);
                end
            end
        end
        diff=max(max(abs(mHelp-mValueFunction1)));
        mValueFunction1 = mHelp;        
    end	
    %mPolicyAsset = vGridAsset(mPolicyAsset);
end
end