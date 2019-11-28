function [mAsset,mConsumption] = Simulation_InfiniteHorizon(mPolicyAsset,...
    vGridAsset,vGridShock,nHousehold,nPeriod)


nGridShock = length(vGridShock);

rng(1)
shockDraws = normrnd(0,sigmaY,[nPeriod,nHousehold]);

mIncome = zeros(nPeriod,nHousehold);
mLogIncome = zeros(nPeriod,nHousehold);
mLogIncome(:,1) = vGridLogShock((nGridShock+1)/2);

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
    
    mAsset(t,hh) = mPolicyAsset(indAsset,indIncome);
    mConsumption(t-1,hh) = mIncome(t-1,hh) + (1+r)*mAsset(t-1,hh)...
        - mAsset(t,hh);
end
end