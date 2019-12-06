%% Housekeeping

close all
outpath = '/Users/Leon/Desktop/PhD-Studies/2nd_Year/Courses/Macro_HetHH/Project/Tex/Figures/Part1_PE/';
%outpath ='Output/';
addpath('Functions')


%% Parameters for Income Process and Grids

logShockAverage = 1.5;
truncOpt = 0;       % Trunctation of normal distribution shock in AR(1) of income process    

ssigma = 1;      % persistence of income shock
ddelta = 0.8;
rrho = 0.04;
r= 0.02;

nGridAsset = 5000;
nGridShock = 31;

minGridAsset = 0;
maxGridAsset = 75;      


%% 4. Infinite Horizon


% Low Variance
ssigmaY = 0.2;
nGridAsset = 5000;
maxGridAsset = 250; 
[vGridAssetLow,vGridShockLow,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

optAccelerator = 10; % Accelerator for VFI, only maximize every two iterations
[it,mValueFunctionInfLow,mPolicyAssetInfLow,mPolicyConsInfLow] = VFI_InfHorizon(rrho,r,ssigma,vGridAssetLow,vGridShockLow,mTransitionShock,0,optAccelerator);


% High Variance
ssigmaY = 0.4;
maxGridAsset = 750; 
[vGridAssetHigh,vGridShockHigh,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

optAccelerator = 10;
[ithigh,mValueFunctionInfHigh,mPolicyAssetInfHigh,mPolicyConsInfHigh] = VFI_InfHorizon(rrho,r,ssigma,vGridAssetHigh,vGridShockHigh,mTransitionShock,0,optAccelerator);


% Plots Infinite Horizon

figure(1)
plot(vGridAssetLow,mPolicyConsInfLow(:,1),vGridAssetLow,mPolicyConsInfLow(:,(nGridShock+1)/2),vGridAssetLow,mPolicyConsInfLow(:,end),'Linewidth',1);
xlabel('assets');
ylabel('consumption');
legend('z_{low}','z_{avg}','z_{high}');
set(gca,'FontSize',13,'Fontweight','bold');
print('-dpng', [outpath,'consFunc_inf_low','.png']);


figure(2)
plot(vGridAssetHigh,mPolicyConsInfHigh(:,1),vGridAssetHigh,mPolicyConsInfHigh(:,(nGridShock+1)/2),vGridAssetHigh,mPolicyConsInfHigh(:,end),'Linewidth',1);
xlabel('assets');
ylabel('consumption');
legend('z_{low}','z_{avg}','z_{high}');
set(gca,'FontSize',13,'Fontweight','bold');
print('-dpng', [outpath,'consFunc_inf_high','.png']);

%{
mConsShareLow = mPolicyConsInfLow./(mPolicyConsInfLow+mPolicyAssetInfLow);
mConsShareHigh = mPolicyConsInfHigh./(mPolicyConsInfHigh+mPolicyAssetInfHigh);

figure(3)
plot(vGridAssetLow,mConsShareLow(:,1)-mConsShareHigh(:,1),'Linewidth',1.0);
xlabel('assets');
ylabel('consumption');
legend('z_{low }');
set(gca,'FontSize',13,'Fontweight','bold');
%print('-depsc', [outpath,'consFunc_inf_high','.eps']);
print('-dpng', [outpath,'consFunc_inf_comp1','.png']);

figure(4)
plot(vGridAssetLow,mConsShareLow(:,end)-mConsShareHigh(:,end),'Linewidth',1.0);
xlabel('assets');
ylabel('consumption');
legend('z_{high}');
set(gca,'FontSize',13,'Fontweight','bold');
%print('-depsc', [outpath,'consFunc_inf_high','.eps']);
print('-dpng', [outpath,'consFunc_inf_comp2','.png']);


figure(5)
plot(vGridAssetLow,mPolicyAssetInfLow(:,1),vGridAssetLow,mPolicyAssetInfLow(:,(nGridShock+1)/2),vGridAssetLow,mPolicyAssetInfLow(:,end),'Linewidth',1.5);
xlabel('assets');
ylabel('savings');
legend('z_{low}','z_{avg}','z_{high}');
set(gca,'FontSize',13,'Fontweight','bold');
%print('-depsc', [outpath,'assetFunc_inf_low','.eps']);
print('-dpng', [outpath,'assetFunc_inf_low','.png']);


figure(6)
plot(vGridAssetHigh,mPolicyAssetInfHigh(:,1),vGridAssetHigh,mPolicyAssetInfHigh(:,(nGridShock+1)/2),vGridAssetHigh,mPolicyAssetInfHigh(:,end),'Linewidth',1.5);
xlabel('assets');
ylabel('savings');
legend('z_{low}','z_{avg}','z_{high}');
set(gca,'FontSize',13,'Fontweight','bold');
%print('-depsc', [outpath,'assetFunc_inf_high','.eps']);
print('-dpng', [outpath,'assetFunc_inf_high','.png']);
%}

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


%% 5. Finite Horizon

nGridAsset = 5000;
nPeriod = 61;
MortOpt = 0;


% Low Variance
ssigmaY = 0.2;
maxGridAsset = 75; 
[vGridAsset,vGridShockLow,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

[mValueFunctionFiniteLow,mPolicyAssetFiniteLow,mPolicyConsFiniteLow] = VFI_FinHorizon(rrho,r,ssigma,vGridAsset,vGridShockLow,mTransitionShock,nPeriod,MortOpt);


% High Variance
ssigmaY = 0.4;
maxGridAsset = 75; 
[vGridAsset,vGridShockHigh,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

[mValueFunctionFiniteHigh,mPolicyAssetFiniteHigh,mPolicyConsFiniteHigh] = VFI_FinHorizon(rrho,r,ssigma,vGridAsset,vGridShockHigh,mTransitionShock,nPeriod,MortOpt);


% Plots Finite Horizon

% Consumption function for young and old agents with medium income
figure(3)
plot(vGridAsset,mPolicyConsFiniteLow(:,(nGridShock+1)/2,end),vGridAsset,mPolicyConsFiniteLow(:,(nGridShock+1)/2,25),'Linewidth',1.5);
xlabel('assets');
ylabel('consumption');
legend('young','old');
set(gca,'FontSize',13,'Fontweight','bold');
%print('-dpng', [outpath,'consFunc_fin_low','.png']);


figure(4)
plot(vGridAsset,mPolicyConsFiniteHigh(:,(nGridShock+1)/2,end),vGridAsset,mPolicyConsFiniteHigh(:,(nGridShock+1)/2,25),'Linewidth',1.5);
xlabel('assets');
ylabel('consumption');
legend('young','old');
set(gca,'FontSize',13,'Fontweight','bold');
%print('-dpng', [outpath,'consFunc_fin_high','.png']);

%{
figure(3)
plot(vGridAsset,F1(vGridAsset',6*ones(1,nGridAsset)'),vGridAsset,F2(vGridAsset',6*ones(1,nGridAsset)'),'Linewidth',1.5);
xlabel('assets');
ylabel('consumption');
legend('z_{low}','z_{high}');
set(gca,'FontSize',13,'Fontweight','bold');
%print('-depsc', [outpath,'consFunc_inf_low','.eps']);
print('-dpng', [outpath,'consFunc_fin_comp','.png']);


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



%% 6.1 Simulation for 60 Periods without borrowing

ssigma = 1;
ssigmaY = 0.2;
nGridAsset = 5000;
minGridAsset = 0;
maxGridAsset = 75; 

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);


nPeriod = 61;
MortOpt = 0;

[~,mPolicyAsset,mPolicyCons] = VFI_FinHorizon(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


nHousehold = 1e04;
IncomeDataOpt = 0;

[mAsset,mConsumption,~] = Simulation_FiniteHorizon(mPolicyAsset,mPolicyCons,vGridAsset,vGridShock,ssigmaY,ddelta,nHousehold,IncomeDataOpt);


vAvgConsPath = mean(mConsumption,2);
vAvgAssetPath = mean(mAsset,2);
time = 20:80;

figure;
pl=plot(time,vAvgConsPath,time,vAvgAssetPath);
xla=xlabel('Age');
yla=ylabel('Asset and Consumption');
%tit=title('Average Consumption and Asset over lifecycle without Data');
le=legend('Consumption','Asset');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
%set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'FontSize',14,'Fontweight','bold');
print('-dpng', [outpath,'sim_simple','.png']);



% For differrent risk aversion
ssigma = 3;
[~,mPolicyAsset,mPolicyCons] = VFI_FinHorizon(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);
[mAsset,mConsumption,~] = Simulation_FiniteHorizon(mPolicyAsset,mPolicyCons,vGridAsset,vGridShock,ssigmaY,ddelta,nHousehold,IncomeDataOpt);
vAvgConsPathNew = mean(mConsumption,2);
vAvgAssetPathNew = mean(mAsset,2);

time = 20:79;
figure;
pl=plot(time,vAvgConsPath(1:60),time,vAvgConsPathNew(1:60));
xla=xlabel('Age');
yla=ylabel('Consumption');
%tit=title('Average Consumption over lifecycle with different risk aversion');
le=legend('ssigma=1','ssigma=3');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
%set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'FontSize',14,'Fontweight','bold');
print('-dpng', [outpath,'sim_risk','.png']);






%% 6.2 Simulation for 60 Periods with borrowing 

% High Variance
ssigmaY = 0.2;
ssigma = 1;
nGridAsset = 350;
minGridAsset = -2;
maxGridAsset = 75; 

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);


nPeriod = 61;
MortOpt = 0;

[mValueFunction,mPolicyAsset,mPolicyCons] = VFI_FinHorizon_Bor(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


nHousehold = 1e04;
IncomeDataOpt = 0;

[mAsset,mConsumption,~] = Simulation_FiniteHorizon_Bor(mPolicyAsset,mPolicyCons,vGridAsset,vGridShock,ssigmaY,ddelta,nHousehold,IncomeDataOpt);


vAvgConsPath = mean(mConsumption,2);
vAvgAssetPath = mean(mAsset,2);
time = 20:80;

figure;
pl=plot(time,vAvgConsPath,time,vAvgAssetPath);
xla=xlabel('Age');
yla=ylabel('Asset and Consumption');
%tit=title('Average Consumption and Asset over lifecycle without Data');
le=legend('Consumption','Asset');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
%set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'FontSize',14,'Fontweight','bold');
print('-dpng', [outpath,'sim_borrowing','.png']);




%% 7. Simulation for 60 Periods with Income and Mortality Data

% High Variance
ssigmaY = 0.2;
ssigma = 1;
nGridAsset = 250;
minGridAsset = -15;
maxGridAsset = 75;      

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);
 
nPeriod = 61;
MortOpt = 1;

[~,mPolicyAsset,mPolicyCons] = VFI_FinHorizon_Bor(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


nHousehold = 1e04;
IncomeDataOpt = 1;

[mAsset,mConsumption,~] = Simulation_FiniteHorizon_Bor(mPolicyAsset,mPolicyCons,vGridAsset,vGridShock,ssigmaY,ddelta,nHousehold,IncomeDataOpt);


vAvgConsPath = mean(mConsumption,2);
vAvgAssetPath = mean(mAsset,2);
time = 20:80;

figure;
pl=plot(time,vAvgConsPath,time,vAvgAssetPath);
xla=xlabel('Age');
yla=ylabel('Asset and Consumption');
%tit=title('Average Consumption and Assets over lifecycle with Data');
le=legend('Consumption','Asset');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
%set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'FontSize',14,'Fontweight','bold');
print('-dpng', [outpath,'sim_Data','.png']);




%% 8.1 Comparison with Empirical Life Cycle Consumption Profile

% Model 
ssigmaY = 0.2;
ssigma = 1;
nGridAsset = 10000;
maxGridAsset = 100;  
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);


nPeriod = 61;
MortOpt = 1;

[~,mPolicyAsset,mPolicyCons] = VFI_FinHorizon(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


nHousehold = 1e04;
IncomeDataOpt = 1;

[mAsset,mConsumption,mIncome] = Simulation_FiniteHorizon(mPolicyAsset,mPolicyCons,vGridAsset,vGridShock,ssigmaY,ddelta,nHousehold,IncomeDataOpt);


vAvgConsPath = mean(mConsumption,2);
%vAvgAssetPath = mean(mAsset,2);
time = 20:79;

% Data
inpath = 'Data/';
mConsData = load([inpath,'consprofile.txt']); 

vConsData = sum(reshape(mConsData(:,2),[4,length(mConsData)/4]),1);
vConsData = vConsData*(vAvgConsPath(1)/vConsData(1));
vAge = mConsData(1,1):mConsData(end,1);


figure;
pl=plot(time,vAvgConsPath(1:60),time,vConsData(1:60));
xla=xlabel('Age');
yla=ylabel('Consumption');
%tit=title('Average Consumption over lifecycle Model vs Data');
le=legend('Model','Data');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
%set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'FontSize',14,'Fontweight','bold');
print('-dpng', [outpath,'Data_comparison','.png']);




%% 8.2 Comparison with Empirical Life Cycle Consumption Profile

% Model 
ssigmaY = 0.03;
ssigma = 1;
nGridAsset = 10000;
maxGridAsset = 100;  
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);


nPeriod = 61;
MortOpt = 1;

[~,mPolicyAsset,mPolicyCons] = VFI_FinHorizon(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,nPeriod,MortOpt);


nHousehold = 1e04;
IncomeDataOpt = 1;

[mAsset,mConsumption,mIncome] = Simulation_FiniteHorizon(mPolicyAsset,mPolicyCons,vGridAsset,vGridShock,ssigmaY,ddelta,nHousehold,IncomeDataOpt);


vAvgConsPath = mean(mConsumption,2);
%vAvgAssetPath = mean(mAsset,2);
time = 20:79;

% Data
inpath = 'Data/';
mConsData = load([inpath,'consprofile.txt']); 

vConsData = sum(reshape(mConsData(:,2),[4,length(mConsData)/4]),1);
vConsData = vConsData*(vAvgConsPath(1)/vConsData(1));
vAge = mConsData(1,1):mConsData(end,1);


figure;
pl=plot(time,vAvgConsPath(1:60),time,vConsData(1:60));
xla=xlabel('Age');
yla=ylabel('Consumption');
%tit=title('Average Consumption over lifecycle Model vs Data');
le=legend('Model','Data');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
%set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'FontSize',14,'Fontweight','bold');
set(le,'FontSize',14,'Fontweight','bold');
print('-dpng', [outpath,'Data_comparison2','.png']);





%% 9. Consumption Insurance

ssigmaY = 0.2;
ssigma = 1;
nGridAsset = 1500;
nPeriod = 2000;

nHousehold = 1e02;
optAccelerator = 10;

ddelta = 0.98;
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);
[~,~,mPolicyAsset,~] = VFI_InfHorizon(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,0,optAccelerator);
[~,mConsumption,mIncome] = Simulation_InfiniteHorizon(r,ddelta,ssigmaY,mPolicyAsset,vGridAsset,vGridShock,nHousehold,nPeriod);

mConsumption = mConsumption(1000:end-1,:);
mIncome = mIncome(1000:end-1,:);
insCoeff = zeros(1,nHousehold);

for nh=1:nHousehold
    diffLogC = diff(log(mConsumption(:,nh)));
    diffLogY = diff(log(mIncome(:,nh)));
    covs = cov(diffLogC,diffLogY);
    covs = covs(1,2);
    insCoeff(nh) = 1-covs/var(diffLogY);
end
insCoeff_High = mean(insCoeff);
disp(insCoeff_High)


ddelta = 0.01;
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);
[~,~,mPolicyAsset,~] = VFI_InfHorizon(rrho,r,ssigma,vGridAsset,vGridShock,mTransitionShock,0,optAccelerator);
[~,mConsumption,mIncome] = Simulation_InfiniteHorizon(r,ddelta,ssigmaY,mPolicyAsset,vGridAsset,vGridShock,nHousehold,nPeriod);

mConsumption = mConsumption(1000:end-1,:);
mIncome = mIncome(1000:end-1,:);
insCoeff = zeros(1,nHousehold);

for nh=1:nHousehold
    diffLogC = diff(log(mConsumption(:,nh)));
    diffLogY = diff(log(mIncome(:,nh)));
    covs = cov(diffLogC,diffLogY);
    covs = covs(1,2);
    insCoeff(nh) = 1-covs/var(diffLogY);
end
insCoeff_Low = mean(insCoeff);
disp(insCoeff_Low)




%%

toc;
beep;