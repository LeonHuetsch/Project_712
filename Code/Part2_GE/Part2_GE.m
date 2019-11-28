tic;

addpath('Functions')
outpath ='Output/';

% Choose productivity A such that deterministic SS capital = 1, then
% generate Grid for asset holdings a around SS capital.


%% Aiyagari's Parameters for Income Process and Grids

bbeta = 0.96;
rrho = 1/bbeta - 1;
aalpha = 0.36;
depreciation = 0.08;
A = 1;

%vSsigma = [1,3,5];
vDdelta = [0.6,0.6,0.9];  % Persistence of income shock
%vSsigmaY = [0.2,0.4];   % Variance of income shock    

nGridAsset = 250;
nGridShock = 21;

minGridAsset = 0;
maxGridAsset = 10;

logShockAverage = 0;
truncOpt = 0;

%bDifference=zeros(length(vDdelta),length(vSsigma),length(vSsigmaY));
%bInterestRate=zeros(length(vDdelta),length(vSsigma),length(vSsigmaY));
%bExitFlag=zeros(length(vDdelta),length(vSsigma),length(vSsigmaY));


%% GE tryout
tic
for d=1:length(vDdelta)
ddelta=vDdelta(d);    
vSsigma = [1,3,5];
vSsigmaY = [0.2,0.4];

    for sy=1:length(vSsigmaY)
    ssigmaY=vSsigmaY(sy);
    
        for s=1:length(vSsigma)
        ssigma=vSsigma(s);
        
            capitalDemand = @(r) (aalpha*A/(r+depreciation))^(1/(1-aalpha));
            wage = @(r) (1-aalpha)*A*(aalpha*A/(r+depreciation))^(aalpha/(1-aalpha));


            % Plot capital vs asset holdings (demand and supply)
            %{
            nInterest = 10;
            optAccelerator = 10;

            %vGridInterest = linspace(-depreciation+1e-08,rho-1e-08,nInterest);
            %vGridInterest = linspace(-depreciation+0.15,0.039,nInterest);
            vGridInterest = linspace(-0.01,0.038,nInterest);
            vCapitalDemand = zeros(nInterest,1);
            vExpectedAssetNext = zeros(nInterest,1);

            [vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
                nGridAsset,minGridAsset,maxGridAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

            for i=1:nInterest
                [~,mValueFunction,mPolicyAsset,~] = VFI_InfHorizon(rrho,vGridInterest(i),ssigma,...
                    aalpha,A,depreciation,vGridAsset,vGridShock,mTransitionShock,0,optAccelerator);

                [~,vExpectedAssetNext(i)] = StationaryDist(vGridAsset,nGridShock,...
                mTransitionShock,mPolicyAsset);
                vCapitalDemand(i) = capitalDemand(vGridInterest(i));
            end
            %}
            %{
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
            %}

            [~,indexGuess] = min(abs(vCapitalDemand-vExpectedAssetNext));

            tol=1e-01;
            difference=100;
            it=1;

            lb=-depreciation+0.001;
            ub=rrho-0.001;

            options = optimset('fzero');
            options.TolX = 1e-05;

            findr = @(r) SavingsGivenR(vGridAsset,vGridShock,mTransitionShock,nGridShock,rrho,r,aalpha,A,depreciation,ssigma,mValueFunction)-capitalDemand(r);
            [interestRateRCE,diff,exitflag] = fzero(findr,0.03,options);




            % Calculate Stationary Distribution in RCE
            %{
            [vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
                nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt);

            [mValueFunction,mPolicyAsset,mPolicyCons,~] = ...
                Infinite_horizon_iteration(rrho,interestRateRCE,aalpha,A,depreciation,sigma,...
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
            %}

            %bDifference(d,s,sy)=diff;
            %bInterestRate(d,s,sy)=interestRateRCE*100;
            %bExitFlag(d,s,sy)=exitflag;

        end
    end
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