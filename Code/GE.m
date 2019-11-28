%% Parameters

tic;

outpath ='Output/';
addpath('VFI_Matrix')

% Choose productivity A such that deterministic SS capital = 1, then
% generate Grid for asset holdings a around SS capital.
r = 0.01208;
Kss = 1;
depreciation = 1;
alpha = 0.66;
%A = (r+depreciation)/alpha;
A = 1;
CapDemand = (alpha*A/(r+depreciation))^(1/(1-alpha));

sigmaY = 0.4;   % variance of income shock
delta = 0.8;    % persistence of income shock

nGridAsset = 100;
nGridShock = 21;

logShockAverage = 0;

sigma = 1;      % risk aversion
rho = 0.04;     % discount rate
beta = 1/(1+rho);


%% Grids for z and k, transition matrix for z

m = (nGridShock-1)/2;
logShockMin = logShockAverage - m*sqrt(sigmaY^2/(1-delta^2));
logShockMax = logShockAverage + m*sqrt(sigmaY^2/(1-delta^2));
vGridLogShock=linspace(logShockMin,logShockMax,nGridShock);

% Caldulate transition matrix for shocks
mTransitionShock = zeros(nGridShock);
for j=1:nGridShock
for i=1:nGridShock
    mTransitionShock(i,j) = normcdf((vGridLogShock(j)+0.5*sigmaY/(sqrt(1-delta^2))...
                -(1-delta)*logShockAverage-delta*vGridLogShock(i))/sigmaY)...
             -normcdf((vGridLogShock(j)-0.5*sigmaY/(sqrt(1-delta^2))...
                -(1-delta)*logShockAverage-delta*vGridLogShock(i))/sigmaY);
end
end

for j=1:nGridShock
        mTransitionShock(j,1) = normcdf((vGridLogShock(1)+0.5*sigmaY/(sqrt(1-delta^2))...
                    -(1-delta)*logShockAverage-delta*vGridLogShock(j))/sigmaY);
        mTransitionShock(j,nGridShock) = 1 - normcdf((vGridLogShock(nGridShock)-0.5*sigmaY/(sqrt(1-delta^2))...
                    -(1-delta)*logShockAverage-delta*vGridLogShock(j))/sigmaY);
end

vGridShock=exp(vGridLogShock); % the ln() follows an AR(1), so need to exp

%vGridAsset = exp(linspace(log(0.5*Kss),log(100*Kss),nAsset));
vGridAsset = linspace(0,30*Kss,nGridAsset);




%% Value Function Iteration

bAssetToday = repmat(reshape(vGridAsset,[nGridAsset,1]),[1,nGridShock,nGridAsset]);
bAssetNext = repmat(reshape(vGridAsset,[1,1,nGridAsset]),[nGridAsset,nGridShock,1]);
bShockToday = repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1,nGridAsset]);

bConsumption = bShockToday + (1+r)*bAssetToday - bAssetNext;
bConsumption(bConsumption<0) = 0;

maxit = 1e04;
tol = 1e-08;
diff = 100;
it = 0;

mValueFunction = zeros(nGridAsset,nGridShock);

while it<=maxit && diff>tol
    bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
        [1,nGridShock,nGridAsset]),[nGridAsset,1,1]);
    bValue = log(bConsumption) + beta*bContinuation;
    bValue(bConsumption<0) = -1e20;
    
    mHelp=max(bValue,[],3);
    
    diff=max(max(abs(mHelp-mValueFunction)));    
    mValueFunction=mHelp;   
    it=it+1;
end		 

[~,mIndexPolicyAsset]=max(bValue,[],3);
mPolicyAsset=vGridAsset(mIndexPolicyAsset);
mPolicyCons=bConsumption(:,:,mIndexPolicyAsset);


%% Some Plots 
%{
figure;

pl=mesh(mValueFunction);
xla=ylabel('Assets Today');
yla=xlabel('Shock Today');
tit=title('ValueFunction');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');


figure;

pl=mesh(mPolicyAsset);
xla=ylabel('Assets Today');
yla=xlabel('Shock Today');
tit=title('Policy Function');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');


figure;

pl=plot(vGridAsset,mValueFunction(:,1));
xla=ylabel('Assets Today');
yla=xlabel('Shock Today');
tit=title('Policy Function');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');


figure;

pl=plot(vGridAsset,mValueFunction(:,(nGridShock+1)/2));
xla=ylabel('Assets Today');
yla=xlabel('Shock Today');
tit=title('Policy Function');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');


figure;

pl=plot(vGridAsset,mValueFunction(:,end));
xla=ylabel('Assets Today');
yla=xlabel('Shock Today');
tit=title('Policy Function');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');

%}

%toc;



%{
%% Simulation

nHousehold = 1000;    % Number of households
nPeriod = 61;      % Number of Periods

rng(1)
shockDraws = normrnd(0,sigmaY,[nPeriod,nHousehold]);

mIncome = zeros(nPeriod,nHousehold);
mLogIncome = zeros(nPeriod,nHousehold);
mLogIncome(1,:) = vGridLogShock((nGridShock+1)/2);
mIncome(1,:) = exp(mLogIncome(1,:));

mConsumption = zeros(nPeriod,nHousehold);
mAsset = zeros(nPeriod,nHousehold);
mAsset(1,:) = vGridAsset(nGridAsset/2);
%mAsset(1,:) = vGridAsset;

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

figure
mesh(mAsset)

figure
mesh(mConsumption)
%}




%% Stationary Distribution

% Transition matrix (a1y1,a1y2,..., a2y1,a1y2,...) for measure Phi
nStates = nGridShock*nGridAsset;
mTransitionState=zeros(nStates);
for yc=1:nGridShock
for ac=1:nGridAsset
        ind=find(vGridAsset==mPolicyAsset(ac,yc));
        mTransitionState(yc+(ac-1)*nGridShock,(ind-1)*nGridShock...
            +1:(ind-1)*nGridShock+nGridShock) = mTransitionShock(yc,:);
end
end

% First eigenvector has eigenvalue 1, iteration is faster though
%[mEigenvector,vEigenvalue] = eig(mTransitionState');
%vEigenvector1 = mEigenvector(:,1);
%vStationaryDist = vEigenvector1/sum(vEigenvector1);

vStationaryDist = ones(nStates,1)/(nStates);
for i=1:1000
    vStationaryDist = mTransitionState' * vStationaryDist;
end

mStationaryDist = reshape(vStationaryDist,[nGridShock,nGridAsset]);
% Expected asset holdings next period
expectedAssetNext = sum(sum(mStationaryDist'.*mPolicyAsset));

%mesh(mStationaryDist)
disp(CapDemand)
disp(expectedAssetNext)

% Put all of this in a function to minimize




%% Check Uniqueness

capitalDemand = @(r) (alpha*A/(r+depreciation))^(1/(1-alpha));

nInterest = 50;
%vGridInterest = linspace(-depreciation+1e-08,rho-1e-08,nInterest);
vGridInterest = linspace(-0.6,0.025,nInterest);
vCapital = zeros(nInterest,1);
vExpectedAssetNext = zeros(nInterest,1);

for i=1:nInterest
    [~,mPolicyAsset,~,~] = Infinite_horizon_iteration(rho,vGridInterest(i),...
        vGridAsset,vGridShock,mTransitionShock);
    mStationaryDist = StationaryDist(vGridAsset,nGridShock,...
    mTransitionShock,mPolicyAsset);

    vExpectedAssetNext(i) = sum(sum(mStationaryDist'.*mPolicyAsset));
    vCapital(i) = capitalDemand(vGridInterest(i));
end


figure;

pl=plot(vGridInterest,vCapital,vGridInterest,vExpectedAssetNext);
xla=ylabel('Interest Rate');
tit=title('Capital Stock and Expected Asset Holdings');
le=legend('Capital Stock','Expected Asset Holdings');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(le,'Fontsize',12,'Fontweight','bold');
print('-depsc', [outpath,'Capital_AssetHoldings','.eps']);

toc;