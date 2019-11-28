tic;
%bInterestRates = zeros(3,2,4);
%bDifferences = zeros(3,2,4);

vdelta = [0,0.3,0.6,0.9];

for d=1:length(vdelta)
    
    vsigmaY = [0.2,0.4];
    vsigma = [1,3,5];
    delta=vdelta(d);
    vHelpDiff = zeros(3,2,1);
    vHelpInterest = zeros(3,2,1);    
    
    for sy=1:length(vsigmaY)
        sigmaY=vsigmaY(sy);
    for s=1:length(vsigma)
        sigma=vsigma(s);
        
        beta = 0.96;
        rho = 1/beta - 1;
        alpha = 0.36;
        depreciation = 0.08;
        A = 1;
        kappa = 0.5;
        
        nGridAsset = 50;
        nGridShock = 15;
        nGridLabor = 2;

        minGridAsset = 0;
        maxGridAsset = 25;

        logShockAverage = 0;

        ub = rho;
        lb = -depreciation;
        diff1=100;
        tol1=1e-02;
        maxit1=50;
        it=1;
        
        
        while diff1>tol1 && it<=maxit1
            
            it=it+1;
            r = (ub+lb)/2;
        
            %% Grid and Transition Matrix for Income Process    
            vHelp=linspace(-3,3,nGridShock);
            vGridLogShock = logShockAverage + sqrt(sigmaY^2/(1-delta^2))*vHelp;
            stepSize=(vGridLogShock(nGridShock)-vGridLogShock(1))/(nGridShock-1);

            mTransitionShock = zeros(nGridShock);    
            for j=1:nGridShock
            for i=1:nGridShock
                if j == 1
                    mTransitionShock(i,j) = normcdf((vGridLogShock(1) - delta*vGridLogShock(i)...
                                            + stepSize/2) / sigmaY);
                elseif j == nGridShock
                    mTransitionShock(i,j) = 1 - normcdf((vGridLogShock(nGridShock) - ...
                                            delta*vGridLogShock(i) - stepSize/2) / sigmaY);
                else
                    mTransitionShock(i,j) = normcdf((vGridLogShock(j) - delta*vGridLogShock(i)...
                                            + stepSize/2) / sigmaY) - ...
                                            normcdf((vGridLogShock(j) - delta*vGridLogShock(i)...
                                            - stepSize/2) / sigmaY);
                end
            end
            end

            vGridShock=exp(vGridLogShock); % the ln() follows an AR(1), so need to exp

            % Normalize income process, so average income is 1
            [mEigenVector,mEigenValue]=eig(mTransitionShock');
            [~,indexEV]=min(abs(diag(real(mEigenValue))-1));
            vStaionaryDist=mEigenVector(:,indexEV)/sum(mEigenVector(:,indexEV));
            vGridShock = vGridShock/(vGridShock*vStaionaryDist);


            %vGridAsset = exp(linspace(log(0.01),log(maxGridAsset),nGridAsset));
            vGridAsset = linspace(minGridAsset,maxGridAsset,nGridAsset);



            %% Value Function Iteration
            wage = (1-alpha)*A*(alpha*A/(r+depreciation))^(alpha/(1-alpha));

            beta = 1/(1+rho);
            nGridAsset = length(vGridAsset);
            nGridShock = length(vGridShock);

            bAssetToday = repmat(reshape(vGridAsset,[nGridAsset,1]),...
                [1,nGridShock,nGridAsset*nGridLabor]);
            bShockToday = repmat(reshape(vGridShock,[1,nGridShock]),...
                [nGridAsset,1,nGridAsset*nGridLabor]);
            bAssetNext = repmat(reshape(vGridAsset,[1,1,nGridAsset]),...
                [nGridAsset,nGridShock,nGridLabor]);
            bLaborToday = zeros(nGridAsset,nGridShock,nGridAsset*nGridLabor);
            bLaborToday(:,:,1:nGridAsset) = 1;
          
            bConsumption = wage*bShockToday.*bLaborToday +...
                (1+r)*bAssetToday - bAssetNext;

            bConsumption(bConsumption<=0) = 0;

            maxit2 = 1e04;
            tol2 = 1e-5;
            diff2 = 100;
            it2 = 0;

            %if mValueGuess == 0
            mValueFunction = zeros(nGridAsset,nGridShock);
            %else
                %mValueFunction = mValueGuess;
            %end

            if sigma == 1
                bUtility = log(bConsumption) - kappa*bLaborToday;
                while it2<=maxit2 && diff2>tol2
                    bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
                        [1,nGridShock,nGridAsset]),[nGridAsset,1,nGridLabor]);
                    bValue = bUtility + beta*bContinuation;
                    bValue(bConsumption<=0) = -1e20;

                    mHelp=max(bValue,[],3);

                    diff2=max(max(abs(mHelp-mValueFunction)));    
                    mValueFunction=mHelp;   
                    it2=it2+1;
                end	

            else
                bUtility = ((bConsumption.^(1-sigma))-1)/(1-sigma) - kappa*bLaborToday;
                while it2<=maxit2 && diff2>tol2
                    bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
                        [1,nGridShock,nGridAsset]),[nGridAsset,1,1]);
                    bValue = bUtility + beta*bContinuation;
                    bValue(bConsumption<0) = -1e20;

                    mHelp=max(bValue,[],3);

                    diff2=max(max(abs(mHelp-mValueFunction)));    
                    mValueFunction=mHelp;   
                    it2=it2+1;
                end	    
            end

            [~,mIndexPolicy]=max(bValue,[],3);
            
            mIndexPolicyAsset = rem(mIndexPolicy,nGridAsset);
            mIndexPolicyAsset(mIndexPolicyAsset==0)=nGridAsset;
            mPolicyAsset=vGridAsset(mIndexPolicyAsset);
            
            mPolicyLabor=mIndexPolicy<=nGridAsset;
            %mPolicyCons = wage*bShockToday(:,:,1) + (1+r)*bAssetToday(:,:,1)...
            %        - mPolicyAsset;


            %% Stationary distribution
            nStates = nGridShock*nGridAsset;

            % Transition matrix (a1y1,a1y2,..., a2y1,a1y2,...) for measure Phi
            mTransitionState=zeros(nStates);
            for ac=1:nGridAsset
                for yc=1:nGridShock
                    ind=find(vGridAsset==mPolicyAsset(ac,yc));
                    mTransitionState(yc+(ac-1)*nGridShock,(ind-1)*nGridShock...
                        +1:(ind-1)*nGridShock+nGridShock) = mTransitionShock(yc,:);
                end
            end

            % First eigenvector has eigenvalue 1, iteration is faster though
            %[mEigenvector,vEigenvalue] = eig(mTransitionState');
            %vEigenvector1 = mEigenvector(:,1);
            %vStationaryDist = vEigenvector1/sum(vEigenvector1);

            vStationaryDist = ones(nStates,1)/(nStates);
            maxit3=1e04;
            it3=1;
            tol3=1e-06;
            diff3=100;

            while it3<maxit3 && diff3>tol3
                vHelp = mTransitionState' * vStationaryDist;

                diff3 = max(abs(vHelp-vStationaryDist));
                vStationaryDist = vHelp;
                it3=it3+1;
            end

            mStationaryDist = reshape(vStationaryDist,[nGridShock,nGridAsset]);


            %% Equilibrium 

            expectAssetHoldings = sum(sum(mStationaryDist'.*mPolicyAsset));
            capitalDemand = (alpha*A/(r+depreciation))^(1/(1-alpha));

            difference = capitalDemand-expectAssetHoldings;
            
            if difference<0
                ub = r;
            else
                lb = r;
            end
            diff1=abs(difference);
                
        end
        vHelpDiff(s,sy)=diff1;
        vHelpInterest(s,sy)=r;
    end
    end
    %bInterestRates(:,:,d) = vHelpInterest;
    %bDifferences(:,:,d) = vHelpDiff;
    %toc;
end
toc;