 function [vGridAsset,vGridShock,mTransitionShock] = SetupGrids(...
    nGridAsset,minGridAsset,maxGridAsset,nGridShock,sigmaY,delta,logShockAverage,truncOpt)

% Function generates Grids for assets and income shocks as well as the
% transition matrix for the income process.

% For truncOpt = 0, the innovation term of the AR(1) is not
% truncated, normal with 0 mean, sigmaY^2 variance
% For truncOpt != 0, the innovation term of the AR(1) follows
% a normal with 0 mean, sigmaY^2 variance and symmetrically
% truncated at the truncOpt value


if truncOpt == 0
    
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
    
    
else
    
    truncValue = abs(truncOpt);
    vHelp=linspace(-3,3,nGridShock);
    vGridLogShock = logShockAverage + sigmaY/sqrt(1-delta^2)*vHelp;
    
    
    % Caldulate transition matrix for shocks
    mTransitionShock = zeros(nGridShock);
    
    vGridLogShockExtend = [-Inf,vGridLogShock,Inf];
    cdf_trunc = @(x,mean,std) max(min((cdf('Normal',x,mean,std)...
        - cdf('Normal',mean-truncValue*std,mean,std))/(cdf('Normal',mean+truncValue*std,mean,std)...
        - cdf('Normal',mean-truncValue*std,mean,std)),1),0);
    
    for j=1:nGridShock
        mTransitionShock(j,:) = cdf_trunc((vGridLogShockExtend(2:end-1)+vGridLogShockExtend(3:end))/2,...
            delta*vGridLogShock(j)+(1-delta)*logShockAverage,sigmaY)...
            - cdf_trunc((vGridLogShockExtend(1:end-2)+vGridLogShockExtend(2:end-1))/2,...
            delta*vGridLogShock(j)+(1-delta)*logShockAverage,sigmaY);
    end

    vGridShock=exp(vGridLogShock); % the ln() follows an AR(1), so need to exp
    
    % Normalize income process, so average income is 1
    [mEigenVector,d]=eig(mTransitionShock');
    [~,indexEV]=min(abs(diag(real(d))-1));
    vStaionaryDist=mEigenVector(:,indexEV)/sum(mEigenVector(:,indexEV));
    vGridShock = vGridShock/(vGridShock*vStaionaryDist);

    
    %vGridAsset = exp(linspace(log(minGridAsset),log(maxGridAsset),nGridAsset));
    vGridAsset = linspace(minGridAsset,maxGridAsset,nGridAsset);    
end
end
