dBins=[[25, 100, 300];[100, 300, inf]];
predCorrThreshs = 0:0.1:0.5;
overFitThreshs = 0.1:.05:.2;
nShuffles = 0;
rawBeta = nan(length(predCorrThreshs),length(overFitThreshs),...
    size(dBins,2),3,3); %predCorr by overFit by dBin by inf-type by 3-preds
normBeta = rawBeta; %normalized by std of design matrix
sigBeta = rawBeta; %bootstrap significance of coefficient

noCorrInd = isnan(sigCorr) | isnan(noiCorr);
opt = glmnetSet(struct('alpha',1e-2,'nlambda',250,'lambda_min',1e-5));
for iPred = 1:length(predCorrThreshs)
    iPred,
    predCorrThresh = predCorrThreshs(iPred);
    for iFit = 1:length(overFitThreshs)
        iFit,
        overFitThresh = overFitThreshs(iFit);
        overFitInd = (nTRC-nTC)>overFitThresh | (tTRC-tTC)>overFitThresh;
        badFitInd = nTC<predCorrThresh | tTC<predCorrThresh;
        for dBin = 1:size(dBins,2)
            dBin,
            tmpInd = iD>dBins(1,dBin) & iD<dBins(2,dBin) & ~overFitInd & ~badFitInd & ~noCorrInd;
            X = [iD(tmpInd),sigCorr(tmpInd),noiCorr(tmpInd)];
            Yn = sN(tmpInd);
            Yp = sP(tmpInd);
            Yr = rM(tmpInd)./sqrt(rV(tmpInd));
            
            % remove Nan values from regression-based influence
            Xr = X(~isnan(Yr),:);
            Yr = Yr(~isnan(Yr),:);
            
            % Yn
%             fprintf('\n Yn Shuffle'),
            bShuf = nan(size(X,2),nShuffles);
            parfor nShuffle = 1:nShuffles
                shufOrder = randperm(size(X,1));
                data = cvglmnet(X(shufOrder,:),Yn,'gaussian',opt,'deviance',20);
                bShuf(:,nShuffle) = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
            end
%             fprintf('\n Yn Fit'),
            data = cvglmnet(X,Yn,'gaussian',opt,'deviance',20);
            bOpt = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
            rawBeta(iPred,iFit,dBin,1,:) = bOpt;
            normBeta(iPred,iFit,dBin,1,:) = bOpt./nanstd(X)';
            sigBeta(iPred,iFit,dBin,1,:) = mean(bOpt > bShuf,2);
            
            % Yp
%             fprintf('\n Yp Shuffle'),
            bShuf = nan(size(X,2),nShuffles);
            parfor nShuffle = 1:nShuffles
                shufOrder = randperm(size(X,1));
                data = cvglmnet(X(shufOrder,:),Yp,'gaussian',opt,'deviance',20);
                bShuf(:,nShuffle) = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
            end
%             fprintf('\n Yp Fit'),
            data = cvglmnet(X,Yp,'gaussian',opt,'deviance',20);
            bOpt = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
            rawBeta(iPred,iFit,dBin,2,:) = bOpt;
            normBeta(iPred,iFit,dBin,2,:) = bOpt./nanstd(X)';
            sigBeta(iPred,iFit,dBin,2,:) = mean(bOpt > bShuf,2);
            
            % Yr
%             fprintf('\n Yr Shuffle'),
            bShuf = nan(size(Xr,2),nShuffles);
            parfor nShuffle = 1:nShuffles
                shufOrder = randperm(size(Xr,1));
                data = cvglmnet(Xr(shufOrder,:),Yr,'gaussian',opt,'deviance',20);
                bShuf(:,nShuffle) = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
            end
%             fprintf('\n Yr Fit'),
            data = cvglmnet(Xr,Yr,'gaussian',opt,'deviance',20);
            bOpt = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
            rawBeta(iPred,iFit,dBin,3,:) = bOpt;
            normBeta(iPred,iFit,dBin,3,:) = bOpt./nanstd(Xr)';
            sigBeta(iPred,iFit,dBin,3,:) = mean(bOpt > bShuf,2);
        end
    end
end

% save('170922_corrByDistanceInfluence','rawBeta','normBeta','sigBeta','predCorrThreshs','overFitThreshs','dBins'),

%%
thisMat = rawBeta;
% thisMat = -log10(sigBeta); thisMat(isinf(thisMat))=max(thisMat(isfinite(thisMat)));

allNear = reshape(squeeze(thisMat(:,:,1,:,:)),[],3);
allMid = reshape(squeeze(thisMat(:,:,2,:,:)),[],3);
allFar = reshape(squeeze(thisMat(:,:,3,:,:)),[],3);
X = [allNear(:,2:3),allMid(:,2:3),allFar(:,2:3)];

figure,imagesc(X),
[e,s,l]=pca(X);

figure,plot(e(1:2:end,1)),hold on,plot(e(2:2:end,1)),
title(sprintf('First Eig: %2.1f%s total var',100*l(1)/sum(l),'%')),
xticks(1:3), xticklabels({'25-100um','100-300um','>300um'})
legend('Signal Correlation','Noise Correlation')

figure,imagesc(reshape(s(:,1),length(predCorrThreshs),[])),
title(sprintf('First Eig: %2.1f%s total var',100*l(1)/sum(l),'%')),
yticks(1:6),yticklabels({'0','0.1','0.2','0.3','0.4','0.5'}),