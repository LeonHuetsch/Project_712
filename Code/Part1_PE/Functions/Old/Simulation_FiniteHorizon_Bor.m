function [mAsset,mConsumption,mIncome] = Simulation_FiniteHorizon_Bor(mPolicyAsset,mPolicyCons,vGridAsset,vGridShock,sigmaY,delta,nHousehold,IncomeDataOpt)

% Simulates consumption asset holding paths for nHousehold households and
% nPeriod periods. All households start with average income and no assets
% If DataOpt=1 it includes income data from incprofile.txt

minAsset = vGridAsset(1);
maxAsset = vGridAsset(end);

if IncomeDataOpt == 0
    
    % Make 3rd dimension go from young to old
    mPolicyCons = mPolicyCons(:,:,end:-1:1);
    mPolicyAsset = mPolicyAsset(:,:,end:-1:1);

    nPeriod = size(mPolicyCons,3);
    nGridShock = length(vGridShock);
    vGridLogShock = log(vGridShock);

    rng(1)
    shockDraws = normrnd(0,sigmaY,[nPeriod,nHousehold]);

    mIncome = zeros(nPeriod,nHousehold);                % Income matrix
    mLogIncome = zeros(nPeriod,nHousehold);             % Log income matrix
    %mLogIncome(1,:) = vGridLogShock((nGridShock+1)/2);  % Start of with medium income
    mLogIncome(1,:) = vGridLogShock(1);
    mIncome(1,:) = exp(mLogIncome(1,:));

    mConsumption = zeros(nPeriod,nHousehold);
    [~,indz] = min(abs(vGridAsset));
    mAsset = zeros(nPeriod,nHousehold);
    mAsset(1,:) = vGridAsset(indz);                   % start with zero assets

    for hh=1:nHousehold    
        for t=2:nPeriod
            
            % Find next income shock on Grid
            [~,IndexNext] = min(abs(vGridLogShock - (delta*mLogIncome(t-1,hh)...
                + sqrt(1-delta^2)*shockDraws(t-1,hh))));
            mLogIncome(t,hh) = vGridLogShock(IndexNext);
            mIncome(t,hh) = exp(mLogIncome(t,hh));

            indIncome = find(vGridLogShock==mLogIncome(t-1,hh));

            
            mAsset(t,hh) = interp1(vGridAsset',mPolicyAsset(:,indIncome,t-1),mAsset(t-1,hh));

            mConsumption(t-1,hh) = interp1(vGridAsset',mPolicyCons(:,indIncome,t-1),mAsset(t-1,hh));

        end
    end

    
else
    
    theta = 0.7;
    
    inpath = 'Data/';
    vIncomeData = load([inpath,'incprofile.txt']);    
    vIncomeData(46:61) = theta*vIncomeData(45);
    
    %[~,indRetire] = min(abs(vGridShock - theta*vIncomeData(45)));
    
    % Make 3rd dimension go from young to old
    mPolicyCons = mPolicyCons(:,:,end:-1:1);
    mPolicyAsset = mPolicyAsset(:,:,end:-1:1);

    nPeriod = size(mPolicyCons,3);
    nGridShock = length(vGridShock);
    vGridLogShock = log(vGridShock);

    rng(1)
    shockDraws = normrnd(0,sigmaY,[nPeriod,nHousehold]);

    mIncome = zeros(nPeriod,nHousehold);                % Income matrix
    mLogIncome = zeros(nPeriod,nHousehold);             % Log income matrix
    mLogIncome(1,:) = vGridLogShock((nGridShock+1)/2);  % Start of with medium income
    %mLogIncome(1,:) = vGridLogShock(1);
    mIncome(1,:) = exp(mLogIncome(1,:));

    mConsumption = zeros(nPeriod,nHousehold);
    [~,indz] = min(abs(vGridAsset));
    mAsset = zeros(nPeriod,nHousehold);
    mAsset(1,:) = vGridAsset(indz+1);  

    for hh=1:nHousehold    
    for t=2:nPeriod
        
        %if t<=45
            % Find next income shock on Grid
            [~,IndexNext] = min(abs(vGridLogShock - (log(vIncomeData(t)) +...
                delta*mLogIncome(t-1,hh) + sqrt(1-delta^2)*shockDraws(t-1,hh))));
            mLogIncome(t,hh) = vGridLogShock(IndexNext);
            mIncome(t,hh) = exp(mLogIncome(t,hh));
        %else
        %    IndexNext = indRetire;
        %    mLogIncome(t,hh) = vGridLogShock(IndexNext);
        %    mIncome(t,hh) = exp(mLogIncome(t,hh));
        %end
        
        indIncome = find(vGridLogShock==mLogIncome(t-1,hh));
        mAsset(t,hh) = min(maxAsset,interp1(vGridAsset',mPolicyAsset(:,indIncome,t-1),mAsset(t-1,hh)));

        mConsumption(t-1,hh) = interp1(vGridAsset',mPolicyCons(:,indIncome,t-1),mAsset(t-1,hh));
    end
    end
end

end