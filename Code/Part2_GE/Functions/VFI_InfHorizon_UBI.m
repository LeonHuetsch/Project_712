function [it,mValueFunction1,mPolicyAsset1,mPolicyCons1,mPolicyLabor1] = ...
    VFI_InfHorizon_UBI(tau,kappa,lambda,rrho,r,ssigma,aalpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,mValueGuess,optAccelerator)


wage = (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));
nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);
bbeta = 1/(1+rrho);

maxit = 1e06;
tol = 1e-05;
diff = 100;
it = 0;

if mValueGuess == 0
    mValueFunction1w = zeros(nGridAsset,nGridShock);
    mValueFunction1s = zeros(nGridAsset,nGridShock);
else
    mValueFunction1 = mValueGuess;
end
mPolicyAsset1w = zeros(nGridAsset,nGridShock);
mPolicyCons1w = zeros(nGridAsset,nGridShock);
mHelpw = zeros(nGridAsset,nGridShock);

mPolicyAsset1s = zeros(nGridAsset,nGridShock);
mPolicyCons1s = zeros(nGridAsset,nGridShock);
mHelps = zeros(nGridAsset,nGridShock);

if ssigma == 1
    while it<=maxit && diff>tol
        it=it+1;
        mContValue = mValueFunction1w*mTransitionShock';

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
                                 

                        %workVal = (log((1-tau)*wage*shockToday+(1+r)*assetToday+lambda-assetNext) - kappa);
                        %leisureVal = log((1+r)*assetToday+lambda-assetNext);
                        cons = (1-tau)*wage*shockToday + (1+r)*assetToday + lambda - assetNext;

                        %if cons<=0.0 
                        %    value=-1e06;
                        %else
                            value=(log(cons)-kappa) + bbeta*mContValue(assetNextIndex,shockIndex); 
                        %end
                        if value>valuePrev
                            valuePrev = value;
                            assetChoice = assetNextIndex;
                            assetChoicePrev = assetNextIndex;
                            consChoice = cons;
                            %laborChoice = labor;
                        else      
                            break
                        end
                    end
                    mHelpw(assetTodayIndex,shockIndex) = valuePrev;
                    mPolicyAsset1w(assetTodayIndex,shockIndex) = assetChoice; 
                    mPolicyCons1w(assetTodayIndex,shockIndex) = consChoice;
                    %mPolicyLabor1(assetTodayIndex,shockIndex) = laborChoice;
                else
                    mHelpw(assetTodayIndex,shockIndex) = (1-bbeta)*(log(mPolicyCons1w(assetTodayIndex,shockIndex))-kappa)... %kappa*mPolicyLabor1(assetTodayIndex,shockIndex))...
                        + bbeta*mContValue(mPolicyAsset1w(assetTodayIndex,shockIndex),shockIndex);            
                end
            end
        end
        diff=max(max(abs(mHelpw-mValueFunction1w)));
        mValueFunction1w = mHelpw;        
    end		
    %mPolicyAsset = vGridAsset(mPolicyAsset);
    
    diff = 100;
    it = 0;
    while it<=maxit && diff>tol
        it=it+1;
        mContValue = mValueFunction1s*mTransitionShock';

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
                                 

                        workVal = (log((1-tau)*wage*shockToday+(1+r)*assetToday+lambda-assetNext) - kappa);
                        leisureVal = log((1+r)*assetToday+lambda-assetNext);
                        cons = (1+r)*assetToday + lambda - assetNext;

                        %if cons<=0.0 
                        %    value=-1e06;
                        %else
                            value=log(cons) + bbeta*mContValue(assetNextIndex,shockIndex); 
                        %end
                        if value>valuePrev
                            valuePrev = value;
                            assetChoice = assetNextIndex;
                            assetChoicePrev = assetNextIndex;
                            consChoice = cons;
                            %laborChoice = labor;
                        else      
                            break
                        end
                    end
                    mHelps(assetTodayIndex,shockIndex) = valuePrev;
                    mPolicyAsset1s(assetTodayIndex,shockIndex) = assetChoice; 
                    mPolicyCons1s(assetTodayIndex,shockIndex) = consChoice;
                    %mPolicyLabor1(assetTodayIndex,shockIndex) = laborChoice;
                else
                    mHelps(assetTodayIndex,shockIndex) = (1-bbeta)*(log(mPolicyCons1s(assetTodayIndex,shockIndex)))...
                        + bbeta*mContValue(mPolicyAsset1s(assetTodayIndex,shockIndex),shockIndex);            
                end
            end
        end
        diff=max(max(abs(mHelps-mValueFunction1s)));
        mValueFunction1s = mHelps;        
    end		
    %mPolicyAsset = vGridAsset(mPolicyAsset);
    choice = mValueFunction1w>=mValueFunction1s;
    mValueFunction1 = choice.*mValueFunction1w + (1-choice).*mValueFunction1s;
    mPolicyAsset1 = choice.*mPolicyAsset1w + (1-choice).*mPolicyAsset1s;
    mPolicyCons1 = choice.*mPolicyCons1w + (1-choice).*mPolicyCons1s;
    mPolicyLabor1 = choice;
    
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