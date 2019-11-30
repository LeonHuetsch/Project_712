tic;


lambda = 0.2;

%addpath('VFI_Matrix')
addpath('Part2_GE/Functions')
if lambda == 0
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

nGridAsset = 100;   % Gridsize assets
nGridShock = 21;    % Gridsize income process

minGridAsset = 0;   % smallest asset grid point
maxGridAsset = 25;  % largest asset grid point

logShockAverage = 0; 
truncOpt = 0;



capitalDemand = @(r) (aalpha*A/(r+depreciation))^(1/(1-aalpha));
wage = @(r) (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));




%% Plot Capital vs Asset Holdings for different r

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

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);
tic
[it,mValueFunction1,mPolicyAsset1,mPolicyCons1,mPolicyLabor1] = ...
    VFI_InfHorizon_UBI(tauHelp,kkappa,lambdaHelp,rrho,vGridInterest(5),ssigma,aalpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,0,1);
toc
tic
[mValueFunction,mPolicyAsset,mPolicyCons,mIndexPolicyAsset,mPolicyLabor] = VFiteration_UBI(tauHelp,lambdaHelp,kkappa,rrho,vGridInterest(5),...
        aalpha,A,depreciation,ssigma,vGridAsset,vGridShock,mTransitionShock,0);
toc    
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




%% Find RCE interest rate, disutility value kappa and tax rate tau

EqConditions = @(EqParameters) sum(ConditionsGE(nGridAsset,minGridAsset,maxGridAsset,...
    nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt,rrho,aalpha,A,depreciation,ssigma,...
    mValueFunction,lambda,EqParameters(1),EqParameters(2),EqParameters(3)).^2);



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




%% Calculate Stationary Distribution in RCE

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);
toc
[mValueFunction,mPolicyAsset,mPolicyCons,~,mPolicyLabor] = ...
    VFiteration_UBI(EqTau,lambda,EqKappa,rrho,EqInterestRate,aalpha,A,depreciation,ssigma,...
    vGridAsset,vGridShock,mTransitionShock,0);
toc
[mStationaryDist,expectAssetHoldings] = StationaryDist...
    (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset);
toc



%% Output and Plots

% Table for macroeconomic aggregates
effectiveLaborSup = sum(sum(mStationaryDist.*mPolicyLabor.*...
    repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])));
EqCapital = (aalpha*A/(EqInterestRate+depreciation))^(1/(1-aalpha))*effectiveLaborSup;
EqtotalCons = sum(sum(mPolicyCons.*mStationaryDist));
EqtotalOutput = A*EqCapital^aalpha;
EqWage = wage(EqInterestRate);

VarNames = {'Output','Capital','Consumption','Wage','Interest_Rate'};
TransferPolicies = {'lambda = 0'};

Summary = table(EqtotalOutput,EqCapital,EqtotalCons,EqWage,EqInterestRate);
Summary.Properties.RowNames = TransferPolicies;
Summary.Properties.VariableNames = VarNames;
disp(Summary);


% Table for equilibrium conditions
workingShare = sum(sum(mStationaryDist.*mPolicyLabor));
govtRevenue = EqTau*EqWage*(sum(sum(mStationaryDist.*(mPolicyLabor.*...
    repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])))));

CapitalMarket = table(expectAssetHoldings,EqCapital,EqCapital-expectAssetHoldings);
CapitalMarket.Properties.RowNames = {'CapitalMarket'};
CapitalMarket.Properties.VariableNames = {'Ea','K','difference'};
disp(CapitalMarket);

WorkingShare = table(workingShare,0.8,workingShare-0.8);
WorkingShare.Properties.RowNames = {'WorkingShare'};
WorkingShare.Properties.VariableNames = {'WorkingShare','Target','difference'};
disp(WorkingShare);

UBI = table(govtRevenue,lambda,govtRevenue-lambda);
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



toc;
beep;