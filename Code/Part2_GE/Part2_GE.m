%% Housekeeping

close all
outpath = '/Users/Leon/Desktop/PhD-Studies/2nd_Year/Courses/Macro_HetHH/Project/Tex/Figures/Part2_GE/';
%outpath ='Output/';
addpath('Functions')


%% Aiyagari's Parameters for Income Process and Grids

tic
bbeta = 0.96;
rrho = 1/bbeta - 1;
aalpha = 0.36;
depr = 0.08;
A = 1;

vSsigma = [1,3,5];
vDdelta = [0.0 0.3,0.6,0.9];  % Persistence of income shock
vSsigmaY = [0.2,0.4];   % Variance of income shock    

nGridAsset = 750;
nGridShock = 21;
vMultiSteps = [100,750];


logShockAverage = 0;
truncOpt = 0;



%% GE for specific parameter specification

ssigma = 5;
ddelta = 0.8;
ssigmaY = 0.4;

minAsset = 0;
maxAsset = 250;

% Plot capital vs asset holdings (demand and supply)

capitalDemand = @(r) (aalpha*A/(r+depr))^(1/(1-aalpha));
wage = @(r) (1-aalpha)*A*(aalpha*A/(r+depr))^(aalpha/(1-aalpha));


nInterest = 20;

%vGridInterest = linspace(-depreciation+0.01,rrho-1e-08,nInterest);
%vGridInterest = linspace(-depreciation+0.15,0.039,nInterest);
vGridInterest = linspace(-0.01,0.03,nInterest);
vCapitalDemand = zeros(nInterest,1);
vExpectedAssetNext = zeros(nInterest,1);

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

for i=1:nInterest
    vExpectedAssetNext(i) = SavingsFunction(vGridAsset,vMultiSteps,vGridShock,mTransitionShock,nGridShock,rrho,vGridInterest(i),aalpha,A,depr,ssigma,0);
    vCapitalDemand(i) = capitalDemand(vGridInterest(i));
end


figure(1)
subplot(2,1,1)
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

subplot(2,1,2)
plot(vGridInterest,vCapitalDemand-vExpectedAssetNext,'Linewidth',2);
xla=xlabel('Interest Rate');
tit=title('Difference Capital Stock and Expected Asset Holdings');
ax=gca;
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
print('-dpng', [outpath,'Capital_AssetHoldings','.png']);


[~,initialGuess] = min(abs(vCapitalDemand-vExpectedAssetNext));

findr = @(r) SavingsFunction(vGridAsset,vMultiSteps,vGridShock,mTransitionShock,nGridShock,rrho,r,aalpha,A,depr,ssigma,0)-capitalDemand(r);
[interestRateRCE,diff,exitflag] = fzero(findr,vGridInterest(initialGuess));


% Calculate policy functions and stationary distribution in RCE

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

optAccelerator = 20;
[it,mValueFunction,mPolicyAsset,mPolicyCons] = MultigridVFI_InfHorizon(rrho,interestRateRCE,ssigma,aalpha,A,depr,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,0,optAccelerator);

[mStationaryDist,expectAssetHoldings] = StationaryDist(vGridAsset,nGridShock,mTransitionShock,mPolicyAsset);

capitalRCE = capitalDemand(interestRateRCE);
wageRCE = wage(interestRateRCE);

figure;
pl=mesh(mStationaryDist);
xla=xlabel('asset holdings');
yla=ylabel('income shock');
tit=title('Stationary Distribution in RCE');
ax=gca;
set(pl,'Linewidth',2);
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
set(yla,'Fontsize',14,'Fontweight','bold');
print('-depsc', [outpath,'StationaryDist','.eps']);

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
print('-dpng', [outpath,'Policy_func_assets','.png']);

% Find assets at which policy function crosses 45 degree line
vHelp = find(mPolicyAsset(:,end)./(1:nGridAsset)'<=1);
maxAssetNecessary = vGridAsset(vHelp(1));
disp('max Asset needed:')
disp(maxAssetNecessary)
            
            
            

%% GE for different parameter specifications

% With 750 grid points

% 75 for 1, 0.3, 0.2
% 120 for 1, 0.3, 0.4
% 120 for 1, 0.6, 0.2
% 140 for 1, 0.6, 0.4
% 140 for 1, 0.9, 0.2

% 150 for 3, 0.6, 0.2
% 160 for 3, 0.6, 0.4

% 170 for 5, 0.3, 0.4
% 200 for 5, 0.6, 0.4
% 250 for 5, 0.9, 0.4

tic

minAsset = 0;
maxAsset = 100;

bDifference=zeros(length(vDdelta),length(vSsigma),length(vSsigmaY));
bInterestRate=zeros(length(vDdelta),length(vSsigma),length(vSsigmaY));
bExitFlag=zeros(length(vDdelta),length(vSsigma),length(vSsigmaY));

parfor sy=1:length(vSsigmaY)
    ssigmaY = vSsigmaY(sy);
    vSsigma = [1,3,5];
    vDdelta = [0.0 0.3,0.6,0.9];
    
    mDiffHelp = zeros(length(vDdelta),length(vSsigma));
    mInterestHelp = zeros(length(vDdelta),length(vSsigma));
    mExitFlagHelp = zeros(length(vDdelta),length(vSsigma));
    
    options = optimset('fzero');
    options.TolX = 1e-08;
    
    for s=1:length(vSsigma)
    ssigma=vSsigma(s);
    interestRateRCE = 0.035;
    
        for d=1:length(vDdelta)
        ddelta=vDdelta(d);

            capitalDemand = @(r) (aalpha*A/(r+depr))^(1/(1-aalpha));
            %wage = @(r) (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));
    
            [vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);            

            findr = @(r) SavingsFunction(vGridAsset,vMultiSteps,vGridShock,mTransitionShock,nGridShock,rrho,r,aalpha,A,depr,ssigma,0)-capitalDemand(r);
            [interestRateRCE,diff,exitflag] = fzero(findr,interestRateRCE);%,options);
            
            mDiffHelp(d,s) = diff;
            mInterestHelp(d,s) = interestRateRCE*100;
            mExitFlagHelp(d,s) = exitflag;
        end
    end
    bDifference(:,:,sy)=mDiffHelp;
    bInterestRate(:,:,sy)=mInterestHelp;
    bExitFlag(:,:,sy)=mExitFlagHelp;
end
toc

%% Save Data

vRowNames={'0.0','0.3','0.6','0.9'};
vColumnNames={'One','Three','Five'};


tDifference1 = array2table(bDifference(:,:,1),'RowNames',vRowNames,...
    'VariableNames',vColumnNames);
tDifference2 = array2table(bDifference(:,:,2),'RowNames',vRowNames,...
    'VariableNames',vColumnNames);
InterestRatesTable1 = array2table(bInterestRate(:,:,1),'RowNames',vRowNames,...
    'VariableNames',vColumnNames);
InterestRatesTable2 = array2table(bInterestRate(:,:,2),'RowNames',vRowNames,...
    'VariableNames',vColumnNames);

disp(tDifference1) 
disp(tDifference2)
disp(InterestRatesTable1)
disp(InterestRatesTable2)

writetable(InterestRatesTable1,[outpath,'Aiyagari1.csv']);
writetable(InterestRatesTable2,[outpath,'Aiyagari2.csv']);

dlmwrite([outpath,'bDifference.txt'],bDifference,'delimiter','\t','newline','pc');
dlmwrite([outpath,'bInterestRate.txt'],bInterestRate,'delimiter','\t','newline','unix');
dlmwrite([outpath,'bExitFlag.txt'],bExitFlag,'delimiter','\t','newline','pc');


beep;