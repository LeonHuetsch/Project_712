function [mAsset,mConsumption] = Simulation_FiniteHorizon(mPolicyAsset,...
    mPolicyCons,vGridAsset,vGridShock,sigmaY,delta,nHousehold,IncomeDataOpt)

% Simulates consumption asset holding paths for nHousehold households and
% nPeriod periods. All households start with average income and no assets
% If DataOpt=1 it includes income data from incprofile.txt


if IncomeDataOpt == 0
    % Make 3rd dimension go from young to old
    mPolicyCons = mPolicyCons(:,:,end:-1:1);
    mPolicyAsset = mPolicyAsset(:,:,end:-1:1);

    nPeriod = size(mPolicyCons,3);
    nGridShock = length(vGridShock);
    vGridLogShock = log(vGridShock);

    rng(1)
    shockDraws = normrnd(0,sigmaY,[nPeriod,nHousehold]);

    mIncome = zeros(nPeriod,nHousehold);
    mLogIncome = zeros(nPeriod,nHousehold);
    mLogIncome(1,:) = vGridLogShock((nGridShock+1)/2);
    %mLogIncome(1,:) = vGridLogShock(20);
    mIncome(1,:) = exp(mLogIncome(1,:));

    mConsumption = zeros(nPeriod,nHousehold);
    mAsset = zeros(nPeriod,nHousehold);
    %mAsset(1,:) = vGridAsset(nGridAsset/2);

    for hh=1:nHousehold    
    for t=2:nPeriod
        [~,IndexNext] = min(abs(vGridLogShock - (delta*mLogIncome(t-1,hh)...
            + sqrt(1-delta^2)*shockDraws(t-1,hh))));
        mLogIncome(t,hh) = vGridLogShock(IndexNext);
        mIncome(t,hh) = exp(mLogIncome(t,hh));

        indAsset = find(vGridAsset==mAsset(t-1,hh));
        indIncome = find(vGridShock==mIncome(t-1,hh));
        %indAsset = (mAsset(t-1,hh)==vGridAsset);
        %indIncome = (mIncome(t-1,hh)==vGridShock);

        mAsset(t,hh) = mPolicyAsset(indAsset,indIncome,t-1);
        %mConsumption(t-1,hh) = mIncome(t-1,hh) + (1+r)*mAsset(t-1,hh)...
        %    - mAsset(t,hh);
        mConsumption(t-1,hh) = mPolicyCons(indAsset,indIncome,t-1);
    end
    end

    
else
    
    theta = 0.7;
    
    inpath = 'Data/';
    vIncomeData = load([inpath,'incprofile.txt']);    
    vIncomeData(46:61) = theta*vIncomeData(45);
    
    % Make 3rd dimension go from young to old
    mPolicyCons = mPolicyCons(:,:,end:-1:1);
    mPolicyAsset = mPolicyAsset(:,:,end:-1:1);

    nPeriod = size(mPolicyCons,3);
    nGridShock = length(vGridShock);
    vGridLogShock = log(vGridShock);

    rng(1)
    shockDraws = normrnd(0,sigmaY,[nPeriod,nHousehold]);

    mIncome = zeros(nPeriod,nHousehold);
    mLogIncome = zeros(nPeriod,nHousehold);
    mLogIncome(1,:) = vGridLogShock((nGridShock+1)/2);
    %mLogIncome(1,:) = vGridLogShock(20);
    mIncome(1,:) = exp(mLogIncome(1,:));

    mConsumption = zeros(nPeriod,nHousehold);
    mAsset = zeros(nPeriod,nHousehold);
    %mAsset(1,:) = vGridAsset(nGridAsset/2);

    for hh=1:nHousehold    
    for t=2:nPeriod
        [~,IndexNext] = min(abs(vGridLogShock - (log(vIncomeData(t)) +...
            delta*mLogIncome(t-1,hh) + sqrt(1-delta^2)*shockDraws(t-1,hh))));
        mLogIncome(t,hh) = vGridLogShock(IndexNext);
        mIncome(t,hh) = exp(mLogIncome(t,hh));

        indAsset = find(vGridAsset==mAsset(t-1,hh));
        indIncome = find(vGridLogShock==mLogIncome(t-1,hh));
        %indAsset = (mAsset(t-1,hh)==vGridAsset);
        %indIncome = (mIncome(t-1,hh)==vGridShock);

        mAsset(t,hh) = mPolicyAsset(indAsset,indIncome,t-1);
        %mConsumption(t-1,hh) = mIncome(t-1,hh) + (1+r)*mAsset(t-1,hh)...
        %    - mAsset(t,hh);
        mConsumption(t-1,hh) = mPolicyCons(indAsset,indIncome,t-1);
    end
    end
end

end