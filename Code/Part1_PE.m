tic;

addpath('VFI_Matrix')
addpath('VFI_Loops')
addpath('Functions')
% Choose productivity A such that deterministic SS capital = 1, then
% generate Grid for asset holdings a around SS capital.


%% Parameters for Income Process and Grids

r= 0.02;
rrho = 0.04;

ddelta = 0.65;    % persistence of income shock

nGridAsset = 5000;
nGridShock = 21;

minGridAsset = 0;
maxGridAsset = 150;

logShockAverage = 1.5;
truncOpt = 0;


%% 4. Infinite Horizon

r= 0.02;
rrho = 0.04;
ssigma = 1;

% Low Variance
ssigmaY = 0.2;
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

optAccelerator = 2;
[it,mValueFunctionInfLow,mPolicyAssetInfLow,mPolicyConsInfLow] = ...
    VFI_InfHorizon(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,optAccelerator);


% High Variance
ssigmaY = 0.4;
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

optAccelerator = 5;
[ithigh,mValueFunctionInfHigh,mPolicyAssetInfHigh,mPolicyConsInfHigh] = ...
    VFI_InfHorizon(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,optAccelerator);


%% Plots Infinite Horizon
%{
figure;
hold on
for i = 1:nGridShock
 plot(mPolicyConsInfLow(:,i))
end

figure;
hold on
for i = 1:nGridShock
 plot(mPolicyConsInfHigh(:,i))
end

figure;
hold on
for i = 1:nGridShock
 plot(mValueFunctionInfLow(:,i))
end

figure;
hold on
for i = 1:nGridShock
 plot(mValueFunctionInfHigh(:,i))
end

figure;
hold on
for i = 1:nGridShock
 plot(mPolicyAssetInfLow(:,i))
end

figure;
hold on
for i = 1:nGridShock
 plot(mPolicyAssetInfHigh(:,i))
end
%}


%% 5. Iteration with Finite Horizon
%{
nPeriod = 61;
MortOpt = 0;


% Low Variance
sigmaY = 0.2;
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt);

[mValueFunctionFiniteLow,mPolicyAssetFiniteLow,mPolicyConsFiniteLow] = ...
    Finite_horizon_iteration(rho,r,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


% High Variance
sigmaY = 0.4;
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt);

[mValueFunctionFiniteHigh,mPolicyAssetFiniteHigh,mPolicyConsFiniteHigh] = ...
    Finite_horizon_iteration(rho,r,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


%% Plots Finite Horizon

% Compare Policy Consumption for Young and Old
figure;
pl=mesh(mPolicyConsFiniteLow(:,:,end));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Policy Function Consumption Low Variance at age 0');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');

figure;
pl=mesh(mPolicyConsFiniteLow(:,:,3));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Finite Policy Function Consumption Low Variance at age 55');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');

figure;
pl=mesh(mPolicyConsFiniteHigh(:,:,end));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Policy Function Consumption High Variance at age 0');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');

figure;
pl=mesh(mPolicyConsFiniteHigh(:,:,3));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Policy Function Consumption High Variance at age 55');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');


% Compare Policy Asset Holdings for Young and Old
figure;
pl=mesh(mPolicyAssetFiniteLow(:,:,end));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Policy Function Asset Low Variance at age 0');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');

figure;
pl=mesh(mPolicyAssetFiniteLow(:,:,3));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Policy Function Asset Low Variance at age 55');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');


figure;
pl=mesh(mPolicyAssetFiniteHigh(:,:,end));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Policy Function Asset High Variance at age 0');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');

figure;
pl=mesh(mPolicyAssetFiniteHigh(:,:,3));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Policy Function Asset High Variance at age 55');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');


% Compare Value Funtions for Young and Old
figure;
pl=mesh(mValueFunctionFiniteLow(:,:,end));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Value Function Low Variance at age 0');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');

figure;
pl=mesh(mValueFunctionFiniteLow(:,:,3));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Value Function Low Variance at age 55');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');


figure;
pl=mesh(mValueFunctionFiniteHigh(:,:,end));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Value Function High Variance at age 0');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');

figure;
pl=mesh(mValueFunctionFiniteHigh(:,:,3));
xla=xlabel('Shock Today');
yla=ylabel('Assets Today');
tit=title('Value Function High Variance at age 55');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');

%}


%% 6. Simulation for 60 Periods without Income Data
%{
% High Variance
sigmaY = 0.4;
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt);


nPeriod = 61;
MortOpt = 0;

[~,mPolicyAsset,mPolicyCons] = ...
    Finite_horizon_iteration(rho,r,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


nHousehold = 1e03;
IncomeDataOpt = 0;

[mAsset,mConsumption] = Simulation_FiniteHorizon(mPolicyAsset,...
    mPolicyCons,vGridAsset,vGridShock,sigmaY,delta,nHousehold,IncomeDataOpt);


vAvgConsPath = mean(mConsumption,2);
vAvgAssetPath = mean(mAsset,2);
time = 20:80;

figure;
pl=plot(time,vAvgConsPath,time,vAvgAssetPath);
xla=xlabel('Age');
yla=ylabel('Asset and Consumption');
tit=title('Average Consumption and Asset over lifecycle without Data');
le=legend('Consumption','Asset');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'FontSize',14,'Fontweight','bold');

%}


%% 7. Simulation for 60 Periods with Income Data
%{
% High Variance
sigmaY = 0.4;
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt);


nPeriod = 61;
MortOpt = 0;

[~,mPolicyAsset,mPolicyCons] = ...
    Finite_horizon_iteration(rho,r,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


nHousehold = 1e04;
IncomeDataOpt = 1;

[mAsset,mConsumption] = Simulation_FiniteHorizon(mPolicyAsset,...
    mPolicyCons,vGridAsset,vGridShock,sigmaY,delta,nHousehold,IncomeDataOpt);


vAvgConsPath = mean(mConsumption,2);
vAvgAssetPath = mean(mAsset,2);
time = 20:80;

figure;
pl=plot(time,vAvgConsPath,time,vAvgAssetPath);
xla=xlabel('Age');
yla=ylabel('Asset and Consumption');
tit=title('Average Consumption and Assets over lifecycle with Data');
le=legend('Consumption','Asset');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'FontSize',14,'Fontweight','bold');

%}


%% 8. Comparison with Empirical Life Cycle Consumption Profile
%{
% Model 
sigmaY = 0.4;
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt);


nPeriod = 61;
MortOpt = 1;

[~,mPolicyAsset,mPolicyCons] = ...
    Finite_horizon_iteration(rho,r,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


nHousehold = 1e04;
IncomeDataOpt = 1;

[mAsset,mConsumption] = Simulation_FiniteHorizon(mPolicyAsset,...
    mPolicyCons,vGridAsset,vGridShock,sigmaY,delta,nHousehold,IncomeDataOpt);


vAvgConsPath = mean(mConsumption,2);
vAvgAssetPath = mean(mAsset,2);
time = 20:80;


% Data
inpath = 'Data/';
mConsData = load([inpath,'consprofile.txt']); 

vConsData = sum(reshape(mConsData(:,2),[4,length(mConsData)/4]),1);
vConsData = vConsData*(vAvgConsPath(1)/vConsData(1));
vAge = mConsData(1,1):mConsData(end,1);


figure;
pl=plot(time,vAvgConsPath,vAge,vConsData);
xla=xlabel('Age');
yla=ylabel('Consumption');
tit=title('Average Consumption over lifecycle Model vs Data');
le=legend('Model','Data');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'FontSize',14,'Fontweight','bold');

%}


%%

toc;
beep;