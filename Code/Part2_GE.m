tic;

addpath('VFI_Matrix')
addpath('Functions')
outpath ='Output/Output_Part2/';

% Choose productivity A such that deterministic SS capital = 1, then
% generate Grid for asset holdings a around SS capital.


%% Aiyagari's Parameters for Income Process and Grids

beta = 0.96;
rho = 1/beta - 1;
alpha = 0.36;
depreciation = 0.08;
A = 1;
%A = @(r) (r+depreciation)/alpha;

%vsigma = [1,3,5];
%vdelta = [0.6,0.6,0.9];  % Persistence of income shock
%vsigmaY = [0.2,0.4];   % Variance of income shock    

nGridAsset = 100;
nGridShock = 21;

minGridAsset = 0;
maxGridAsset = 5;

logShockAverage = 0;
truncOpt = 0;

bDifference=zeros(length(vdelta),length(vsigma),length(vsigmaY));
bInterestRate=zeros(length(vdelta),length(vsigma),length(vsigmaY));
bExitFlag=zeros(length(vdelta),length(vsigma),length(vsigmaY));


%% GE tryout
%{
for d=1:length(vdelta)
delta=vdelta(d);    
for sy=1:length(vsigmaY)
sigmaY=vsigmaY(sy);
for s=1:length(vsigma)
sigma=vsigma(s);
%}

        capitalDemand = @(r) (alpha*A/(r+depreciation))^(1/(1-alpha));
        wage = @(r) (1-alpha)*A*(alpha*A/(r+depreciation))^(alpha/(1-alpha));


        %% Plot capital vs asset holdings (demand and supply)

        nInterest = 7;
        %vGridInterest = linspace(-depreciation+1e-08,rho-1e-08,nInterest);
        %vGridInterest = linspace(-depreciation+0.15,0.039,nInterest);
        vGridInterest = linspace(-0.01,0.038,nInterest);
        vCapitalDemand = zeros(nInterest,1);
        vExpectedAssetNext = zeros(nInterest,1);

        [vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
            nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt);

        for i=1:nInterest
            [mValueFunction,mPolicyAsset,~,~] = Infinite_horizon_iteration(rho,vGridInterest(i),...
                alpha,A,depreciation,sigma,vGridAsset,vGridShock,mTransitionShock,0);

            [~,vExpectedAssetNext(i)] = StationaryDist(vGridAsset,nGridShock,...
            mTransitionShock,mPolicyAsset);
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


        [~,indexGuess] = min(abs(vCapitalDemand-vExpectedAssetNext));

        tol=1e-01;
        difference=100;
        it=1;

        lb=-depreciation+0.001;
        ub=rho-0.001;

        [interestRateRCE,diff,exitflag] = fsolve(@(r) SavingsGivenR(nGridAsset,minGridAsset,maxGridAsset,...
            nGridShock,sigmaY,delta,logShockAverage,truncOpt,rho,r,alpha,A,...
            depreciation,sigma,mValueFunction),vGridInterest(indexGuess));



        %% Calculate Stationary Distribution in RCE

        [vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
            nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt);

        [mValueFunction,mPolicyAsset,mPolicyCons,~] = ...
            Infinite_horizon_iteration(rho,interestRateRCE,alpha,A,depreciation,sigma,...
            vGridAsset,vGridShock,mTransitionShock,mValueFunction);

        [mStationaryDist,expectAssetHoldings] = StationaryDist...
            (vGridAsset,nGridShock,mTransitionShock,mPolicyAsset);

        capStockSS = capitalDemand(interestRateRCE);
        wageSS = wage(interestRateRCE);

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
        vHelp = find(mPolicyAsset(:,end)./vGridAsset'<=1);
        maxAssetNecessary = vGridAsset(vHelp(1));
        disp('max Asset needed:')
        disp(maxAssetNecessary)


        %bDifference(d,s,sy)=diff;
        %bInterestRate(d,s,sy)=interestRateRCE*100;
        %bExitFlag(d,s,sy)=exitflag;
%{
end
end
end
%}
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

toc;
beep;