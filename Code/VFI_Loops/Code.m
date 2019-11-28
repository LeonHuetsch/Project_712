%% Inputs of the function / Parameters

tic;

% Choose productivity A such that deterministic SS capital = 1, then
% generate Grid for asset holdings a around SS capital.
r = 0.02;
Kss = 1;
depreciation = 0.8;
alpha = 0.66;
A = (r+depreciation)/alpha;

sigmaY = 0.2;   % variance of income shock
delta = 0.8;    % persistence of income shock

nGridAsset = 20;
nGridShock = 15;

logShockAverage = 0;

sigma = 1;      % risk aversion
rho = 0.04;     % discount rate
beta = 1/(1+rho);


% Grids for z and k, transition matrix for z

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
%vGridAsset = linspace(0.5*Kss,500*Kss,nAsset);
vGridAsset = exp(linspace(log(0.5*Kss),log(100*Kss),nGridAsset));




%% Value Function Iteration

maxit = 1e04;
tol = 1e-08;
diff = 100;
it = 0;

mValueFunction = zeros(nGridAsset,nGridShock);
mPolicyAsset = zeros(nGridAsset,nGridShock);
mHelp = zeros(nGridAsset,nGridShock);

while it<=maxit && diff>tol
        for shockToday=1:nGridShock
            capitalChoice=1;
        for assetToday=1:nGridAsset
            vHelpPrev=-1e12;
            mHelp(assetToday,shockToday)=-1e20;
            for assetNext=capitalChoice:nGridAsset
            %for assetNext=1:nAsset
                cons = vGridShock(shockToday) + (1+r)*vGridAsset(assetToday)...
                    - vGridAsset(assetNext);
                if cons<=0.0 
                    vhelp=-10^12;
                else
                    vhelp=log(cons) +...
                        beta*sum(mValueFunction(assetNext,:).*mTransitionShock(shockToday,:));                        
                end
                if vhelp<vHelpPrev
                    break
                end
                if vhelp>mHelp(assetToday,shockToday) 
                    mHelp(assetToday,shockToday)=vhelp;
                    mPolicyAsset(assetToday,shockToday)=vGridAsset(assetNext); 
                    mPolicyCons(assetToday,shockToday)=cons;
                    capitalChoice=assetNext;
                end
                vHelpPrev=vhelp;
            end
        end
        end
        diff=max(max(abs(mHelp-mValueFunction)));
        mValueFunction = mHelp;        
        it=it+1;
end		 




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



toc;


%% Simulation
nHousehold = 100;    % Number of households
nPeriod = 60;      % Number of Periods

rng(1)
shockDraws = normrnd(0,sigmaY,[nPeriod,nHousehold]);

vLogIncome = zeros(nPeriod,nHousehold);
vLogIncome(:,1) = vGridLogShock((nShock+1)/2);

for hh=1:nHousehold    
for t=2:60
    [~,IndexNext] = min(abs(vGridLogShock - (delta*vLogIncome(t-1,hh)...
        + sqrt(1-delta^2)*shockDraws(t-1,hh))));
    vLogIncome(t,hh) = vGridLogShock(IndexNext);
end
end

vIncome = exp(vLogIncome);


%}
toc;
beep;