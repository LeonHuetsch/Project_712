tic;


llambda = 0.2;

%addpath('VFI_Matrix')
addpath('Part2_GE/Functions')
if llambda == 0
    outpath ='Output_Part2/Output_UBI/NoTransfer';
else
    outpath ='Output_Part2/Output_UBI/Transfer';
end

% Why does everyone sit basically at zero assets. Only up to about asset
% holdings of 4-5 when EqCapital is 5. Change logShockAverage maybe or
% would it be the same as changeing sigmaY?


%% Parameters for Income Process and Grids

bbeta = 0.96;
rrho = 1/bbeta - 1;
aalpha = 0.36;
depreciation = 0.08;
A = 1;
kkappa = 0;
%A = @(r) (r+depreciation)/alpha;

ssigma = 1;
ddelta = 0.8;  % Persistence of income shock
ssigmaY = 0.4;   % Variance of income shock    

nGridAsset = 200;   % Gridsize assets
nGridShock = 21;    % Gridsize income process

minAsset = 0;   % smallest asset grid point
maxAsset = 30;  % largest asset grid point

logShockAverage = 0; 
truncOpt = 0;



capitalDemand = @(r) (aalpha*A/(r+depreciation))^(1/(1-aalpha));
wage = @(r) (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));




%% Plot Capital vs Asset Holdings for different r
%{
nInterest = 5;
%vGridInterest = linspace(-depreciation+1e-08,rho-1e-08,nInterest);
%vGridInterest = linspace(-depreciation+0.15,0.039,nInterest);
vGridInterest = linspace(-0.01,0.038,nInterest);
vCapitalDemand = zeros(nInterest,1);
vExpectedAssetNext = zeros(nInterest,1);
%tauHelp = 0;
lambdaHelp = 0.2;
kkappa = 0.638486722667860;
tauHelp = 0.192680040644605;
%tauHelp = 0;

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

for i=1:nInterest
    [mValueFunction,mPolicyAsset,~,mIndexPolicyAsset,~] = VFiteration_UBI(tauHelp,lambdaHelp,kkappa,rrho,vGridInterest(i),...
        aalpha,A,depreciation,ssigma,vGridAsset,vGridShock,mTransitionShock,0);

    [~,vExpectedAssetNext(i)] = StationaryDist(vGridAsset,nGridShock,mTransitionShock,mIndexPolicyAsset);
    vCapitalDemand(i) = capitalDemand(vGridInterest(i));
end

figure(1)
%subplot(2,1,1)
pl=plot(vGridInterest,vCapitalDemand,vGridInterest,vExpectedAssetNext);
xla=xlabel('Interest Rate');
tit=title('Capital Stock and Expected Asset Holdings');
le=legend('Capital Stock','Expected Asset Holdings');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(le,'Fontsize',12,'Fontweight','bold');
%print('-depsc', [outpath,'Capital_AssetHoldings','.eps']);
%{
subplot(2,1,2)
plot(vGridInterest,vCapitalDemand-vExpectedAssetNext,'Linewidth',2);
xla=xlabel('Interest Rate');
tit=title('Difference Capital Stock and Expected Asset Holdings');
ax=gca;
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
%print('-depsc', [outpath,'diff(r)','.eps']);
%}
%}



%% Stationary equilibrium without UBI: tau=lambda=0, kappa equilibirum object such that working population of 80 percent

tic
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

llambda_noUBI = 0;
ttau_noUBI = 0;
EqConditions_noUBI = @(EqParameters) sum(ConditionsGE(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt,rrho,aalpha,A,depreciation,ssigma,0,llambda_noUBI,EqParameters(1),EqParameters(2),ttau_noUBI).^2);

vInitialGuess = [0.035 0.6];
options = optimset('fminsearch');
options.Display = 'final';
options.TolFun = 1e-06;
options.TolX = 1e-06;

[EqParameters_noUBI,diff_noUBI] = fminsearch(EqConditions_noUBI,vInitialGuess,options);
toc

r_noUBI = EqParameters_noUBI(1);
kkappa_noUBI = EqParameters_noUBI(2);

% Corresponding stationary distribution
[mValueFunction_noUBI,mPolicyAsset_noUBI,mPolicyCons_noUBI,mIndexPolicyAsset_noUBI,mPolicyLabor_noUBI] = ...
    VFiteration_UBI(ttau_noUBI,llambda_noUBI,kkappa_noUBI,rrho,r_noUBI,aalpha,A,depreciation,ssigma,...
    vGridAsset,vGridShock,mTransitionShock,0);
[mStationaryDist_noUBI,expectAssetHoldings_noUBI] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mIndexPolicyAsset_noUBI);


%% Output and Plots for results without UBI

% Table for macroeconomic aggregates
effectiveLaborSup_noUBI = sum(sum(mStationaryDist_noUBI.*mPolicyLabor_noUBI.*repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])));
capital_noUBI = (aalpha*A/(r_noUBI+depreciation))^(1/(1-aalpha))*effectiveLaborSup_noUBI;
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
%print('-depsc', [outpath,'StationaryDist','.eps']);

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
%print('-depsc', [outpath,'Policy_func_assets','.eps']);


%% Stationary distribution for Andrew Young's proposal: lambda=0.2, tau equilibrium object

% Note: kappa is taken from the previous section
tic
kkappa = kkappa_noUBI;
llambda_UBI = 0.2;

EqConditions_UBI = @(EqParameters) sum(ConditionsGE(nGridAsset,minAsset,maxAsset,...
    nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt,rrho,aalpha,A,depreciation,ssigma,...
    0,llambda_UBI,EqParameters(1),kkappa,EqParameters(2)).^2);


vInitialGuess = [0.035 0.1];
options = optimset('fminsearch');
options.Display = 'final';
options.TolFun = 1e-06;
options.TolX = 1e-06;

[EqParameters_UBI,diff_UBI] = fminsearch(EqConditions_UBI,vInitialGuess,options);%,options);
toc

r_UBI = EqParameters_UBI(1);
ttau_UBI = EqParameters_UBI(2);

% Corresponding stationary distribution
[mValueFunction_UBI,mPolicyAsset_UBI,mPolicyCons_UBI,mIndexPolicyAsset_UBI,mPolicyLabor_UBI] = ...
    VFiteration_UBI(ttau_UBI,llambda_UBI,kkappa,rrho,r_UBI,aalpha,A,depreciation,ssigma,...
    vGridAsset,vGridShock,mTransitionShock,0);
[mStationaryDist_UBI,expectAssetHoldings_UBI] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mIndexPolicyAsset_UBI);


%{
%[~,Index] = min(abs(vCapitalDemand-vExpectedAssetNext));
%r0 = vGridInterest(Index);
r0 = 0.037146725946504;
kappa0 = 0.638486722667860;
tau0 = 0.192680040644605;
vInitialGuess = [r0 kappa0 tau0];
lowerBound = [0 0.5 0.1];
upperBound = [0.05 0.7 0.25];

options = optimset('fminsearch');%...
%options = optimoptions('fsolve');
%    ,'UseParallel',true);
%options.Algorithm = 'trust-region';
options.Display = 'final';
options.TolFun = 1e-06;
options.TolX = 1e-06;
%options.UseParallel = true;
[EqParameters,diff] = fminsearch(EqConditions,vInitialGuess,options);
%[EqParameters,diff] = fmincon(EqConditions,vInitialGuess,[],[],[],[],...
%    lowerBound,upperBound,[],options);



%EqParameters=[0.035 0.3 0.1];
EqInterestRate = EqParameters(1);
EqKappa = EqParameters(2);
EqTau = EqParameters(3);
%}



%% Output and Plots for results with UBI

% Table for macroeconomic aggregates
effectiveLaborSup_UBI = sum(sum(mStationaryDist_UBI.*mPolicyLabor_UBI.*...
    repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])));
capital_UBI = (aalpha*A/(r_UBI+depreciation))^(1/(1-aalpha))*effectiveLaborSup_UBI;
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
pl=mesh(mStationaryDist_UBI);
xla=xlabel('income shock');
yla=ylabel('asset holdings');
tit=title('Stationary Distribution with UBI');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'Fontsize',14,'Fontweight','bold');
%print('-depsc', [outpath,'StationaryDist','.eps']);

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
%print('-depsc', [outpath,'Policy_func_assets','.eps']);


%% Calculate Stationary Distribution in RCE
%{
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);
toc
[mValueFunction,mPolicyAsset,mPolicyCons,~,mPolicyLabor] = ...
    VFiteration_UBI(EqTau,llambda,EqKappa,rrho,EqInterestRate,aalpha,A,depreciation,ssigma,...
    vGridAsset,vGridShock,mTransitionShock,0);
toc
[mStationaryDist,expectAssetHoldings] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset);
toc



%% Output and Plots

% Table for macroeconomic aggregates
effectiveLaborSup_UBI = sum(sum(mStationaryDist.*mPolicyLabor.*...
    repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])));
capital_UBI = (aalpha*A/(EqInterestRate+depreciation))^(1/(1-aalpha))*effectiveLaborSup_UBI;
totalCons_UBI = sum(sum(mPolicyCons.*mStationaryDist));
totalOutput_UBI = A*capital_UBI^aalpha;
wage_UBI = wage(EqInterestRate);

VarNames = {'Output','Capital','Consumption','Wage','Interest_Rate'};
TransferPolicy = {'lambda = 0'};

Summary_UBI = table(totalOutput_UBI,capital_UBI,totalCons_UBI,wage_UBI,EqInterestRate);
Summary_UBI.Properties.RowNames = TransferPolicy;
Summary_UBI.Properties.VariableNames = VarNames;
disp(Summary_UBI);


% Table for equilibrium conditions
workingShare_UBI = sum(sum(mStationaryDist.*mPolicyLabor));
govtRevenue_UBI = EqTau*wage_UBI*(sum(sum(mStationaryDist.*(mPolicyLabor.*...
    repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])))));

CapitalMarket_UBI = table(expectAssetHoldings,capital_UBI,capital_UBI-expectAssetHoldings);
CapitalMarket_UBI.Properties.RowNames = {'CapitalMarket'};
CapitalMarket_UBI.Properties.VariableNames = {'Ea','K','difference'};
disp(CapitalMarket_UBI);

WorkingShare_UBI = table(workingShare_UBI,0.8,workingShare_UBI-0.8);
WorkingShare_UBI.Properties.RowNames = {'WorkingShare'};
WorkingShare_UBI.Properties.VariableNames = {'WorkingShare','Target','difference'};
disp(WorkingShare_UBI);

UBI = table(govtRevenue_UBI,llambda,govtRevenue_UBI-llambda);
UBI.Properties.RowNames = {'UBI_Transfer'};
UBI.Properties.VariableNames = {'Govt_Revenue','UBI_Lambda','difference'};
disp(UBI);


figure;
pl=mesh(mStationaryDist);
xla=xlabel('income shock');
yla=ylabel('asset holdings');
tit=title('Stationary Distribution in RCE');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'Fontsize',14,'Fontweight','bold');
%print('-depsc', [outpath,'StationaryDist','.eps']);

figure;
pl=plot(vGridAsset,mPolicyAsset(:,end));
xla=xlabel('asset holdings');
yla=ylabel('savings');
tit=title('Policy Function for highest income realization');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'Fontsize',14,'Fontweight','bold');
%print('-depsc', [outpath,'Policy_func_assets','.eps']);
%}


toc;
beep;