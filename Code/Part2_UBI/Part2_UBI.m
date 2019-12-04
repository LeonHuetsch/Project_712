%% Housekeeping

close all

addpath('Functions')
outpath ='Output/';


%% Parameters for Income Process and Grids
tic 

bbeta = 0.96;           % discount factor
rrho = 1/bbeta - 1;     % discount rate
aalpha = 0.36;          % capital share in Cobb-Douglas production function
depr = 0.08;            % depreciation    
A = 1;                  % productivity parameter in Cobb-Douglas production function
%A = @(r) (r+depreciation)/alpha;

ssigma = 1;             % risk aversion parameter households
ddelta = 0.6;           % persistence of income shock
ssigmaY = 0.2;          % variance of income shock    

nGridAsset = 350;       % gridsize assets
nGridShock = 21;        % gridsize income process

minAsset = 0;           % smallest asset grid point
maxAsset = 30;          % largest asset grid point

logShockAverage = 0;    % average shock in logs
truncOpt = 0;           % truncation of normal distribution for shocks (0 is no truncation)



capitalDemand = @(r,laborSupply) (aalpha*A/(r+depr))^(1/(1-aalpha))*laborSupply;
wage = @(r) (1-aalpha)*A*(aalpha*A/(r+depr))^(aalpha/(1-aalpha));




%% Plot Capital vs Asset Holdings for different r


nInterest = 20;
vGridInterest = linspace(-depr+0.04,rrho-0.0001,nInterest);
%vGridInterest = linspace(-0.02,0.038,nInterest);
vCapitalDemand = zeros(nInterest,1);
vSavings = zeros(nInterest,1);      % Capital supply

llambda = 0.2;
kkappa = 0.6;
ttau = 0.15;            
optAccelerator = 35;    % Accelerator for value function iteration
vMultiSteps = [50 350]; % gridsize steps for multigrid algorithm
mGuessVF = 0;           % initial guess for value function 

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

for rr=1:nInterest
    [mValueFunction,~,~,mIndexPolicyAsset,mPolicyLabor]...
        = MultigridVFI_UBI(ttau,kkappa,llambda,rrho,vGridInterest(rr),ssigma,aalpha,A,depr,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,mGuessVF,optAccelerator);

    [mStationaryDist,vSavings(rr)] = StationaryDist(vGridAsset,nGridShock,mTransitionShock,mIndexPolicyAsset);
    effectLaborSupply = (sum(sum(mStationaryDist.*(mPolicyLabor.*repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])))));

    vCapitalDemand(rr) = capitalDemand(vGridInterest(rr),effectLaborSupply);
end


figure(1)
subplot(2,1,1)
pl=plot(vGridInterest,vCapitalDemand,vGridInterest,vSavings);
xla=xlabel('Interest Rate');
tit=title('Capital Stock and Expected Asset Holdings');
le=legend('Capital Stock','Expected Asset Holdings');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(le,'Fontsize',12,'Fontweight','bold');

subplot(2,1,2)
plot(vGridInterest,vCapitalDemand-vSavings,'Linewidth',2);
xla=xlabel('Interest Rate');
tit=title('Difference Capital Stock and Expected Asset Holdings');
ax=gca;
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
print('-depsc', [outpath,'Capital_MarketClearing','.eps']);







%% Stationary equilibrium without UBI: tau=lambda=0, kappa equilibirum object such that working population of 80 percent

llambda_noUBI = 0;
ttau_noUBI = 0;

vMultiSteps = [60 350];
optAccelerator = 35;    
    
EqConditions_noUBI = @(EqParameters) sum(ConditionsGE(nGridAsset,minAsset,maxAsset,vMultiSteps,nGridShock,ssigmaY,ddelta,...
    logShockAverage,truncOpt,rrho,aalpha,A,depr,ssigma,mGuessVF,llambda_noUBI,EqParameters(1),EqParameters(2),ttau_noUBI,optAccelerator).^2);

vInitialGuess = [0.035 0.6];
options = optimset('fminsearch');
options.Display = 'final';
options.TolFun = 1e-06;
options.TolX = 1e-06;

[EqParameters_noUBI,diff_noUBI] = fminsearch(EqConditions_noUBI,vInitialGuess,options);

r_noUBI = EqParameters_noUBI(1);
kkappa_noUBI = EqParameters_noUBI(2);


% Corresponding stationary distribution
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

[mValueFunction_noUBI,mPolicyAsset_noUBI,mPolicyCons_noUBI,mIndexPolicyAsset_noUBI,mPolicyLabor_noUBI] = ...
    MultigridVFI_UBI(ttau_noUBI,kkappa_noUBI,llambda_noUBI,rrho,r_noUBI,ssigma,aalpha,A,depr,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,mGuessVF,optAccelerator);

[mStationaryDist_noUBI,expectAssetHoldings_noUBI] = StationaryDist(vGridAsset,nGridShock,mTransitionShock,mIndexPolicyAsset_noUBI);




%% Output and Plots for results without UBI

% Table for macroeconomic aggregates
effectiveLaborSup_noUBI = sum(sum(mStationaryDist_noUBI.*mPolicyLabor_noUBI.*repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])));
capital_noUBI = (aalpha*A/(r_noUBI+depr))^(1/(1-aalpha))*effectiveLaborSup_noUBI;
totalCons_noUBI = sum(sum(mPolicyCons_noUBI.*mStationaryDist_noUBI));
totalOutput_noUBI = A*capital_noUBI^aalpha;
wage_noUBI = wage(r_noUBI);

VarNames = {'Output','Capital','Consumption','Wage','Interest_Rate'};
TransferPolicy = {'No_UBI'};

Summary_noUBI = table(totalOutput_noUBI,capital_noUBI,totalCons_noUBI,wage_noUBI,r_noUBI);
Summary_noUBI.Properties.RowNames = TransferPolicy;
Summary_noUBI.Properties.VariableNames = VarNames;
disp(Summary_noUBI);


% Table for equilibrium conditions
workingShare_noUBI = sum(sum(mStationaryDist_noUBI.*mPolicyLabor_noUBI));

CapitalMarket_noUBI = table(expectAssetHoldings_noUBI,capital_noUBI,capital_noUBI-expectAssetHoldings_noUBI);
CapitalMarket_noUBI.Properties.RowNames = {'Capital_Market'};
CapitalMarket_noUBI.Properties.VariableNames = {'Ea','K','difference'};
disp(CapitalMarket_noUBI);

WorkingShare_noUBI = table(workingShare_noUBI,0.8,workingShare_noUBI-0.8);
WorkingShare_noUBI.Properties.RowNames = {'Working_Share'};
WorkingShare_noUBI.Properties.VariableNames = {'Working_Share','Target','difference'};
disp(WorkingShare_noUBI);


figure;
pl=mesh(mStationaryDist_noUBI);
xla=xlabel('income shock');
yla=ylabel('asset holdings');
tit=title('Stationary Distribution without UBI');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'Fontsize',14,'Fontweight','bold');
print('-depsc', [outpath,'StationaryDist_noUBI','.eps']);

figure;
pl=plot(vGridAsset,mPolicyAsset_noUBI(:,end));
xla=xlabel('asset holdings');
yla=ylabel('savings');
tit=title('Policy Function for highest income realization without UBI');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'Fontsize',14,'Fontweight','bold');
print('-depsc', [outpath,'Polfun_assets_noUBI','.eps']);









%% Stationary distribution for Andrew Young's proposal: lambda=0.2, tau equilibrium object

% Note: kappa is taken from the previous section

kkappa = kkappa_noUBI;
llambda_UBI = 0.2;

vMultiSteps = [60 350];
optAccelerator = 35;

EqConditions_UBI = @(EqParameters) sum(ConditionsGE(nGridAsset,minAsset,maxAsset,vMultiSteps,nGridShock,ssigmaY,ddelta,...
    logShockAverage,truncOpt,rrho,aalpha,A,depr,ssigma,mGuessVF,llambda_UBI,EqParameters(1),kkappa,EqParameters(2),optAccelerator).^2);

vInitialGuess = [0.035 0.1];
options = optimset('fminsearch');
options.Display = 'final';
options.TolFun = 1e-06;
options.TolX = 1e-06;

[EqParameters_UBI,diff_UBI] = fminsearch(EqConditions_UBI,vInitialGuess,options);

r_UBI = EqParameters_UBI(1);
ttau_UBI = EqParameters_UBI(2);


% Corresponding stationary distribution
[mValueFunction_UBI,mPolicyAsset_UBI,mPolicyCons_UBI,mIndexPolicyAsset_UBI,mPolicyLabor_UBI] = ...
    MultigridVFI_UBI(ttau_UBI,kkappa,llambda_UBI,rrho,r_UBI,ssigma,aalpha,A,depr,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,mGuessVF,optAccelerator);

[mStationaryDist_UBI,expectAssetHoldings_UBI] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mIndexPolicyAsset_UBI);




%% Output and Plots for results with UBI

% Table for macroeconomic aggregates
effectiveLaborSup_UBI = sum(sum(mStationaryDist_UBI.*mPolicyLabor_UBI.*repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])));
capital_UBI = (aalpha*A/(r_UBI+depr))^(1/(1-aalpha))*effectiveLaborSup_UBI;
totalCons_UBI = sum(sum(mPolicyCons_UBI.*mStationaryDist_UBI));
totalOutput_UBI = A*capital_UBI^aalpha;
wage_UBI = wage(r_UBI);

VarNames = {'Output','Capital','Consumption','Wage','Interest_Rate'};
TransferPolicy = {'UBI'};

Summary_UBI = table(totalOutput_UBI,capital_UBI,totalCons_UBI,wage_UBI,r_UBI);
Summary_UBI.Properties.RowNames = TransferPolicy;
Summary_UBI.Properties.VariableNames = VarNames;
disp(Summary_UBI);


% Table for equilibrium conditions
workingShare_UBI = sum(sum(mStationaryDist_UBI.*mPolicyLabor_UBI));
govtRevenue_UBI = ttau_UBI*wage_UBI*(sum(sum(mStationaryDist_UBI.*(mPolicyLabor_UBI.*...
    repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])))));

CapitalMarket_UBI = table(expectAssetHoldings_UBI,capital_UBI,capital_UBI-expectAssetHoldings_UBI);
CapitalMarket_UBI.Properties.RowNames = {'Capital_Market'};
CapitalMarket_UBI.Properties.VariableNames = {'Ea','K','difference'};
disp(CapitalMarket_UBI);

WorkingShare_UBI = table(workingShare_UBI,0.8,workingShare_UBI-0.8);
WorkingShare_UBI.Properties.RowNames = {'Working_Share'};
WorkingShare_UBI.Properties.VariableNames = {'Working_Share','Target','difference'};
disp(WorkingShare_UBI);

UBI = table(govtRevenue_UBI,llambda_UBI,govtRevenue_UBI-llambda_UBI);
UBI.Properties.RowNames = {'UBI_Transfer'};
UBI.Properties.VariableNames = {'Govt_Revenue','UBI_Lambda','difference'};
disp(UBI);


figure;
pl=mesh(mStationaryDist_UBI(1:75,:));
xla=xlabel('income shock');
yla=ylabel('asset holdings');
tit=title('Stationary Distribution with UBI');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'Fontsize',14,'Fontweight','bold');
print('-depsc', [outpath,'StationaryDist_UBI','.eps']);

figure;
pl=plot(vGridAsset,mPolicyAsset_UBI(:,end));
xla=xlabel('asset holdings');
yla=ylabel('savings');
tit=title('Policy Function for highest income realization with UBI');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'Fontsize',14,'Fontweight','bold');
print('-depsc', [outpath,'Polfun_assets_UBI','.eps']);


% Lorenz Curves
vLorenzX_noUBI = zeros(1,nGridAsset);
vLorenzY_noUBI = zeros(1,nGridAsset);
vLorenzX_UBI = zeros(1,nGridAsset);
vLorenzY_UBI = zeros(1,nGridAsset);
for i=1:nGridAsset
    vLorenzX_noUBI(i) = sum(sum(mStationaryDist_noUBI(1:i,:),2));
    vLorenzY_noUBI(i) = sum(sum(mStationaryDist_noUBI(1:i,:),2).*vGridAsset(1:i)')/expectAssetHoldings_noUBI;
    vLorenzX_UBI(i) = sum(sum(mStationaryDist_UBI(1:i,:),2));
    vLorenzY_UBI(i) = sum(sum(mStationaryDist_UBI(1:i,:),2).*vGridAsset(1:i)')/expectAssetHoldings_UBI;    
end

figure;
plot(linspace(0,1,nGridAsset),linspace(0,1,nGridAsset),vLorenzX_noUBI,vLorenzY_noUBI,vLorenzX_UBI,vLorenzY_UBI,'Linewidth',1.5);
xlabel('share of population from lowest to highest wealth');
ylabel('share of wealth');
title('Lorenz Curves with and without UBI');
legend('45 degree','noUBI','UBI');
set(gca,'FontSize',13,'Fontweight','bold');
print('-depsc', [outpath,'Lorenz_UBI_vs_noUBI','.eps']);


toc;
beep;