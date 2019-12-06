%% Housekeeping

close all
outpath = '/Users/Leon/Desktop/PhD-Studies/2nd_Year/Courses/Macro_HetHH/Project/Tex/Figures/Part3_MIT/';
%outpath ='Output/';
addpath('Functions')


%% Parameters

bbeta = 0.96;
rrho = 1/bbeta - 1;
aalpha = 0.36;
depr = 0.08;

ssigma = 1;
ddelta = 0.8;
ssigmaY = 0.2;

llambda = 0;
ttau = 0;
kkappa = 0;

nGridAsset = 300;
nGridShock = 21;
nStates = nGridAsset*nGridShock;

vMultiSteps = [75,nGridAsset];
optAccelerator = 15;  

minAsset = 0;
maxAsset = 60;

logShockAverage = 0;
truncOpt = 0;




%% Stationary equilibrium for A=1

tic
A = 1;
    
EqCondition = @(r) abs(ConditionsGE(nGridAsset,minAsset,maxAsset,vMultiSteps,nGridShock,ssigmaY,ddelta,...
    logShockAverage,truncOpt,rrho,aalpha,A,depr,ssigma,0,r,optAccelerator));

options = optimset('fminbnd');
options.Display = 'final';
options.TolX = 1e-7;

[rSS,diff] = fminbnd(EqCondition,0.03,0.05,options);


% Corresponding stationary distribution
[vGridAsset,vGridShock,mTransitionShock] = SetupGrids(nGridAsset,minAsset,maxAsset,nGridShock,ssigmaY,ddelta,logShockAverage,truncOpt);

[~,mValueFunctionSS,mIndexPolicyAssetSS,mPolicyConsSS] = MultigridVFI_InfHorizon(rrho,rSS,ssigma,aalpha,A,depr,minAsset,maxAsset,mTransitionShock,vGridShock,vMultiSteps,0,optAccelerator);
mPolicyAssetSS = vGridAsset(mIndexPolicyAssetSS);

[mStationaryDistSS,expectAssetHoldingsSS] = StationaryDist(vGridAsset,nGridShock,mTransitionShock,mIndexPolicyAssetSS);
toc

labsupplySS = 1;



%% Extended Path: Shock hits in t=2 here, t=1 is first period

tolDist = 1e-03;
diffDist = 100;
T = 35;
    
while diffDist>tolDist
    dampener = 0.3;

    vProd = ones(1,T);
    vProd(2:11) = 0.9;

    ExoSeqR = rSS*ones(1,T);
    ExoSeqR(2:end) = linspace(rSS+0.5,rSS,T-1);

    EndoSeqR = zeros(1,T);
    EndoSeqR(1) = rSS;
    EndoSeqR(end) = rSS;
    
    vCapLabRatioExo = ones(1,T)*expectAssetHoldingsSS/labsupplySS;
    vRate = zeros(1,T);

    vCapSupply = zeros(1,T);
    vCapSupply(1) = expectAssetHoldingsSS;
    vCapSupply(2) = expectAssetHoldingsSS;          % Capital supply is determined in t=1
    vCapSupply(end) = expectAssetHoldingsSS;

    bStationaryDist = zeros(nGridAsset,nGridShock,T);
    bStationaryDist(:,:,1) = mStationaryDistSS;
    bStationaryDist(:,:,2) = mStationaryDistSS;     % Distribution is determined in t=1
    bStationaryDist(:,:,end) = mStationaryDistSS;

    mStationaryDist = zeros(nStates,T);
    mStationaryDist(:,1) = reshape(mStationaryDistSS',[nStates,1]);
    mStationaryDist(:,2) = mStationaryDist(:,1);
    mStationaryDist(:,end) = mStationaryDist(:,1);

    bValueFunction = zeros(nGridAsset,nGridShock,T);
    bValueFunction(:,:,1) = mValueFunctionSS;
    bValueFunction(:,:,end) = mValueFunctionSS;

    bPolicyAsset = zeros(nGridAsset,nGridShock,T);
    bPolicyAsset(:,:,1) = mPolicyAssetSS;
    bPolicyAsset(:,:,end) = mPolicyAssetSS;

    bPolicyCons = zeros(nGridAsset,nGridShock,T);
    bPolicyCons(:,:,1) = mPolicyConsSS;
    bPolicyCons(:,:,end) = mPolicyConsSS;

    bPolicyAssetIndex = zeros(nGridAsset,nGridShock,T);
    bPolicyAssetIndex(:,:,1) = mIndexPolicyAssetSS;
    bPolicyAssetIndex(:,:,end) = mIndexPolicyAssetSS;


    % Extended Path method

    maxit = 1e04;
    tol = 1e-04;
    diff = 100;
    it = 0;

    tic
    while it<=maxit && diff>tol
        it = it+1;

        % Backwards for value and policy functions
        for t=1:T-2
            %vRate(T-t) = aalpha*vProd(T-t)*(1/vCapLabRatioExo(T-t))^(1-aalpha)-ddelta;
            [bValueFunction(:,:,T-t),bPolicyAsset(:,:,T-t),~,bPolicyAssetIndex(:,:,T-t)] = ...
                VF_transition(bValueFunction(:,:,T-t+1),rrho,ExoSeqR(T-t),aalpha,vProd(T-t),depr,ssigma,vGridAsset,vGridShock,mTransitionShock);
            %[bValueFunction(:,:,T-t),bPolicyAsset(:,:,T-t),~,~,bPolicyLabor(:,:,T-t)] = ...
            %    VF_transition(bValueFunction(:,:,T-t+1),ttau,llambda,kkappa,rrho,vRate(T-t),aalpha,vProd(T-t),depr,ssigma,vGridAsset,vGridShock,mTransitionShock);
        end

        % Forward for distribution
        for t=2:T-1

            % Transition matrix (a1y1,a1y2,..., a2y1,a1y2,...) for states
            mTransitionState=zeros(nStates);
            for yc=1:nGridShock
                for ac=1:nGridAsset
                    ind=mPolicyAssetIndex(ac,yc,t);
                    mTransitionState(yc+(ac-1)*nGridShock,(ind-1)*nGridShock...
                        +1:(ind-1)*nGridShock+nGridShock) = mTransitionShock(yc,:);
                end
            end
            mStationaryDist(:,t+1) = mTransitionState'*mStationaryDist(:,t);
            bStationaryDist(:,:,t+1) = reshape(mStationaryDist(:,t+1),[nGridShock,nGridAsset,1])';

            % Calculate Asset holdings and labor supply 
            vCapSupply(t) = sum(sum(bStationaryDist(:,:,t-1).*bPolicyAsset(:,:,t-1)));
            %vLabSupply(t) = sum(sum(bStationaryDist(:,:,t).*bPolicyLabor(:,:,t).*repmat(reshape(vGridShock,[1,nGridShock]),[nGridAsset,1])));
            
            EndoSeqR(t) = aalpha*vProd(t)*vCapSupply(t)^(aalpha-1) - depr;
        end
        
        %vCapLabRatioEndo = vCapSupply./vLabSupply;
        %diff = max(abs(vCapLabRatioEndo-vCapLabRatioExo));
        %vCapLabRatioExo = dampener*vCapLabRatioEndo + (1-dampener)*vCapLabRatioExo;
        diff = max(abs(EndoSeqR-ExoSeqR));
        ExoSeqR = dampener*EndoSeqR + (1-dampener)*ExoSeqR;

        if mod(it,10)==0 || it == 1
            fprintf('Iteration = %d, diff = %1.7f\n',it,diff);
        end
    end
    toc    


    % Compare distribution in period T to stationary distribution
    diffDist = max(max(abs(bStationaryDist(:,:,end)-mStationaryDistSS)));
    T=T+2;
end


