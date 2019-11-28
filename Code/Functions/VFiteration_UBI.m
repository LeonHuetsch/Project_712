function [mValueFunction,mPolicyAsset,mPolicyCons,mIndexPolicyAsset,mPolicyLabor] = ...
    VFiteration_UBI(tau,lambda,kappa,rho,r,alpha,A,depreciation,sigma,vGridAsset,...
    vGridShock,mTransitionShock,mValueGuess)


% mValueGuess is the inital guess for the value function. Set to zero if
% there is no initial guess. 
nGridLabor = 2;
wage = (1-alpha)*A*(alpha*A/(r+depreciation))^(alpha/(1-alpha));

beta = 1/(1+rho);
nGridAsset = length(vGridAsset);
nGridShock = length(vGridShock);

minGridAsset = vGridAsset(1);
maxGridAsset = vGridAsset(end);

nGridAssetInitial = 7;
if mValueGuess == 0
    mValueFunction = zeros(nGridAssetInitial,nGridShock);
else
    mValueFunction = mValueGuess;
end

%nMultiGrid = floor(1 + log(nGridAsset/5)/log(2));
nGridAssetIter = 4;

while nGridAssetIter<nGridAsset
    
    nGridAssetIter = min(2*nGridAssetIter - 1,nGridAsset);
        
    F = griddedInterpolant(mValueFunction,'cubic');
    mValueFunction = F(repmat(linspace(1,length(vGridAsset),nGridAssetIter)',...
        [1,nGridShock]),repmat(1:nGridShock,[nGridAssetIter,1]));
    
    vGridAsset = linspace(minGridAsset,maxGridAsset,nGridAssetIter);
    
    %F = griddedInterpolant(mValueFunction);
    %mValueFunction = F(repmat(linspace(1,nGridAsset,nGridAssetIter)',[1,nGridShock]),...
    %    repmat(1:nGridShock,[nGridAssetIter,1]));
    
    bAssetToday = repmat(reshape(vGridAsset,[nGridAssetIter,1]),...
        [1,nGridShock,nGridAssetIter*nGridLabor]);
    bShockToday = repmat(reshape(vGridShock,[1,nGridShock]),...
        [nGridAssetIter,1,nGridAssetIter*nGridLabor]);
    bAssetNext = repmat(reshape(vGridAsset,[1,1,nGridAssetIter]),...
        [nGridAssetIter,nGridShock,nGridLabor]);
    bLaborToday = zeros(nGridAssetIter,nGridShock,nGridAssetIter*nGridLabor);
    bLaborToday(:,:,1:nGridAssetIter) = 1;

    bConsumption = (1-tau)*wage*bShockToday.*bLaborToday + (1+r)*bAssetToday...
        - bAssetNext + lambda;

    bConsumption(bConsumption<=0) = 0;

    maxit = 1e04;
    tol = 1e-5;
    diff = 100;
    it = 0;

    %if mValueGuess == 0
    %    mValueFunction = zeros(nGridAsset,nGridShock);
    %else
    %    mValueFunction = mValueGuess;
    %end

    if sigma == 1
        bUtility = log(bConsumption) - kappa*bLaborToday;
        while it<=maxit && diff>tol
            bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
                [1,nGridShock,nGridAssetIter]),[nGridAssetIter,1,nGridLabor]);
            bValue = bUtility + beta*bContinuation;
            bValue(bConsumption<=0) = -1e20;

            mHelp=max(bValue,[],3);

            diff=max(max(abs(mHelp-mValueFunction)));    
            mValueFunction=mHelp;   
            it=it+1;
        end	

    else
        bUtility = ((bConsumption.^(1-sigma))-1)/(1-sigma) - kappa*bLaborToday;
        while it<=maxit && diff>tol
            bContinuation = repmat(reshape(mTransitionShock*mValueFunction',...
                [1,nGridShock,nGridAssetIter]),[nGridAssetIter,1,nGridLabor]);
            bValue = bUtility + beta*bContinuation;
            bValue(bConsumption<0) = -1e20;

            mHelp=max(bValue,[],3);

            diff=max(max(abs(mHelp-mValueFunction)));    
            mValueFunction=mHelp;   
            it=it+1;
        end	    
    end
end

[~,mIndexPolicy]=max(bValue,[],3);

mIndexPolicyAsset = rem(mIndexPolicy,nGridAssetIter);
mIndexPolicyAsset(mIndexPolicyAsset==0)=nGridAssetIter;
mPolicyAsset=vGridAsset(mIndexPolicyAsset);

mPolicyLabor=mIndexPolicy<=nGridAssetIter;

mPolicyCons = (1-tau)*wage*bShockToday(:,:,1).*mPolicyLabor + (1+r)*bAssetToday(:,:,1)...
        - mPolicyAsset + lambda;

end