function [mAsset,mConsumption,mIncome] = Simulation_InfiniteHorizon(r,ddelta,ssigmaY,mPolicyAsset,vGridAsset,vGridShock,nHousehold,nPeriod)


nGridShock = length(vGridShock);
vGridLogShock = log(vGridShock);

rng(1)
shockDraws = normrnd(0,ssigmaY,[nPeriod,nHousehold]);

mIncome = zeros(nPeriod,nHousehold);
mLogIncome = zeros(nPeriod,nHousehold);
mLogIncome(1,:) = vGridLogShock(1);
mIncome(1,:) = exp(mLogIncome(1,:));

mConsumption = zeros(nPeriod,nHousehold);

[~,indz] = min(abs(vGridAsset));
mAsset = zeros(nPeriod,nHousehold);
mAsset(1,:) = vGridAsset(indz);

for hh=1:nHousehold    
for t=2:nPeriod
    [~,IndexNext] = min(abs(vGridLogShock - (ddelta*mLogIncome(t-1,hh)...
        + sqrt(1-ddelta^2)*shockDraws(t-1,hh))));
    mLogIncome(t,hh) = vGridLogShock(IndexNext);
    mIncome(t,hh) = exp(mLogIncome(t,hh));
    
    indAsset = vGridAsset==mAsset(t-1,hh);
    indIncome = vGridLogShock==mLogIncome(t-1,hh);
    %indAsset = (mAsset(t-1,hh)==vGridAsset);
    %indIncome = (mIncome(t-1,hh)==vGridShock);
    
    mAsset(t,hh) = mPolicyAsset(indAsset,indIncome);
    mConsumption(t-1,hh) = mIncome(t-1,hh) + (1+r)*mAsset(t-1,hh)...
        - mAsset(t,hh);
end
end