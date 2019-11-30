%% Housekeeping

close all
addpath('Functions')
outpath ='Output/';




%% Aiyagari's Parameters for Income Process and Grids

tic
bbeta = 0.96;
rrho = 1/bbeta - 1;
aalpha = 0.36;
depreciation = 0.08;
A = 1;

vSsigma = [1,3,5];
vDdelta = [0.3,0.6,0.9];  % Persistence of income shock
vSsigmaY = [0.2,0.4];   % Variance of income shock    

nGridAsset = 500;
nGridShock = 21;
vMultiSteps = [50,150,500];

minAsset = 0;
maxAsset = 15;

logShockAverage = 0;
truncOpt = 0;




%% GE for specific parameter specification

ssigma = 1;
ddelta = 0.9;
ssigmaY = 0.2;


% Plot capital vs asset holdings (demand and supply)

capitalDemand = @(r) (aalpha*A/(r+depreciation))^(1/(1-aalpha));
wage = @(r) (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));


nInterest = 10;

%vGridInterest = linspace(-depreciation+0.01,rrho-1e-08,nInterest);
%vGridInterest = linspace(-depreciation+0.15,0.039,nInterest);
vGridInterest = linspace(-0.01,0.06,nInterest);
vCapitalDemand = zeros(nInterest,1);
vExpectedAssetNext = zeros(nInterest,1);

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

for i=1:nInterest
    vExpectedAssetNext(i) = SavingsFunction(vGridAsset,vMultiSteps,vGridShock,mTransitionShock,nGridShock,rrho,vGridInterest(i),aalpha,A,depreciation,ssigma,0);
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
%print('-depsc', [outpath,'Capital_AssetHoldings','.eps']);

subplot(2,1,2)
plot(vGridInterest,vCapitalDemand-vExpectedAssetNext,'Linewidth',2);
xla=xlabel('Interest Rate');
tit=title('Difference Capital Stock and Expected Asset Holdings');
ax=gca;
set(ax,'FontSize',14,'Fontweight','bold');
set(tit,'Fontsize',14,'Fontweight','bold');
set(xla,'Fontsize',14,'Fontweight','bold');
%print('-depsc', [outpath,'diff(r)','.eps']);


[~,initialGuess] = min(abs(vCapitalDemand-vExpectedAssetNext));

findr = @(r) SavingsFunction(vGridAsset,vMultiSteps,vGridShock,mTransitionShock,nGridShock,rrho,r,aalpha,A,depreciation,ssigma,0)-capitalDemand(r);
[interestRateRCE,diff,exitflag] = fzero(findr,vGridInterest(initialGuess));


% Calculate policy functions and stationary distribution in RCE

[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

optAccelerator = 10;
[it,mValueFunction,mPolicyAsset,mPolicyCons] = MultigridVFI_InfHorizon(rrho,interestRateRCE,ssigma,aalpha,A,depreciation,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,0,optAccelerator);

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

% Find assets at which policy function crosses 45 degree line
vHelp = find(mPolicyAsset(:,end)./(1:nGridAsset)'<=1);
maxAssetNecessary = vGridAsset(vHelp(1));
disp('max Asset needed:')
disp(maxAssetNecessary)
            
            
            

%% GE for different parameter specifications

tic

bDifference=zeros(length(vDdelta),length(vSsigma),length(vSsigmaY));
bInterestRate=zeros(length(vDdelta),length(vSsigma),length(vSsigmaY));
bExitFlag=zeros(length(vDdelta),length(vSsigma),length(vSsigmaY));

parfor sy=1:length(vSsigmaY)
    ssigmaY = vSsigmaY(sy);
    vSsigma = [1,3,5];
    vDdelta = [0.6,0.75,0.9];
    
    mDiffHelp = zeros(length(vDdelta),length(vSsigma));
    mInterestHelp = zeros(length(vDdelta),length(vSsigma));
    mExitFlagHelp = zeros(length(vDdelta),length(vSsigma));
    
    options = optimset('fzero');
    options.TolX = 1e-05;
    
    for s=1:length(vSsigma)
    ssigma=vSsigma(s);
    interestRateRCE = 0.035;
    
        for d=1:length(vDdelta)
        ddelta=vDdelta(d);

            capitalDemand = @(r) (aalpha*A/(r+depreciation))^(1/(1-aalpha));
            %wage = @(r) (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));
    
            [vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);            

            findr = @(r) SavingsFunction(vGridAsset,vMultiSteps,vGridShock,mTransitionShock,nGridShock,rrho,r,aalpha,A,depreciation,ssigma,0)-capitalDemand(r);
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
%{
%% Save Data

vRowNames={'0.3','0.6','0.9'};
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


dlmwrite('bDifference.txt',bDifference,'delimiter','\t','newline','pc');
dlmwrite('bInterestRate.txt',bInterestRate,'delimiter','\t','newline','unix');
dlmwrite('bExitFlag.txt',bExitFlag,'delimiter','\t','newline','pc');
%}

%toc;
%beep;