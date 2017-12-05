function orchestraStimFitsGP(jobIndex, inputDataFilePath, outputDataFilePath)


load(inputDataFilePath), % should contain allResp structure

nNeurons = 1:size(allResp.Y,2);
kFolds = 20;

ClusterInfo.setWallTime('12:00');
ClusterInfo.setMemUsage('2500');
ClusterInfo.setQueueName('mpi');
parpool(24),

%% Construct fit structures / matrices
X = [cosd(allResp.Dir),sind(allResp.Dir),...
    allResp.SF,allResp.TF,sqrt(allResp.spd)];
Y = sqrt(allResp.Y(:,nNeurons)); Y = Y./repmat(std(Y),size(Y,1),1); 

hyp = struct;
likFunc = {@likGauss};
infFunc = @infExact;
hyp.lik = 0;
meanFunc = [];
covFun1 = {@covMask, {[1 1 0 0 0], {@covSE,'iso',[]}}};
covFun2 = {@covMask, {[0 0 1 0 0], {@covSE,'iso',[]}}}';
covFun3 = {@covMask, {[0 0 0 1 0], {@covSE,'iso',[]}}}';
covFun4 = {@covMask, {[0 0 0 0 1], {@covSE,'iso',[]}}}';
covFunc = {@covScale, {@covProd, {covFun1, covFun2, covFun3, covFun4}}};
hyp.cov = zeros(eval(feval(covFunc{:})),1);

params.meanFunc = meanFunc;
params.covFunc = covFunc;
params.likFunc = likFunc;
params.infFunc = infFunc;
params.initHyp = hyp;
params.X = X;
params.Y = Y;
params.nNeurons = nNeurons;

%% Fit full model
fprintf('\n \n Fitting model to Full dataset \n'),

trainPreds = nan(size(Y));
trainProbs = nan(size(Y));
optHyp = cell(size(Y,2));
parfor nNeuron = 1:size(Y,2)
    if sum(Y(:,nNeuron)) > 1
        optHyp{nNeuron} = minimize(hyp, @gp, -100, ...
            infFunc, meanFunc, covFunc, likFunc, X, Y(:,nNeuron));
        [trainPreds(:,nNeuron),~,~,~,trainProbs(:,nNeuron)] = gp(optHyp{nNeuron},infFunc,...
            meanFunc,covFunc,likFunc, X, Y(:,nNeuron), X, Y(:,nNeuron));
    end
end

optHyp = cat(1,optHyp{:});
if length(optHyp) ~= size(Y,2)
    warning('Some Cells were ignored due to no response data'),
    validFits = find(sum(Y)>1);
    tmpHyp(validFits) = optHyp;
    optHyp = tmpHyp;
end

%% Fit cross-validated models

cvFolds = cvpartition(size(Y,1),'KFold',kFolds);
allPreds = nan(size(Y));
allProbs = nan(size(Y));
for k = 1:kFolds
    fprintf('\n \n Predicting fold %d of %d\n',k,kFolds),
    xTrain = X(cvFolds.training(k),:);
    yTrain = Y(cvFolds.training(k),:);
    xTest = X(cvFolds.test(k),:);
    yTest = Y(cvFolds.test(k),:);
    foldPreds = nan(size(yTest));
    foldProbs = nan(size(yTest));
    parfor nNeuron = 1:size(Y,2)
        if sum(Y(:,nNeuron)) > 1
            [foldPreds(:,nNeuron),~,~,~,foldProbs(:,nNeuron)] = gp(optHyp(nNeuron), infFunc, meanFunc, ...
                covFunc, likFunc, xTrain, yTrain(:,nNeuron), xTest, yTest(:,nNeuron));
        end
    end
    allPreds(cvFolds.test(k),:) = foldPreds;      
    allProbs(cvFolds.test(k),:) = foldProbs; 
end

params.testPreds = allPreds;
for nNeuron = 1:size(Y,2)
    params.predCorr(nNeuron,:) = corr(Y(:,nNeuron),...
        [allPreds(:,nNeuron),trainPreds(:,nNeuron)]);
end

params.trainProbs = trainProbs;
params.testProbs = allProbs;
%% Optimize contrast-phase parameters

xPh = [cosd(3:3:360)',sind(3:3:360)'];
yPh = zscore(allResp.yPh(:,nNeurons));

hypPh = struct;
covFuncPh = @covSEiso; hypPh.cov = [0; 0];
likFuncPh = @likGauss; hypPh.lik = 0;

params.covFuncPh = covFuncPh;
params.likFuncPh = likFuncPh;
params.initHypPh = hypPh;
params.xPh = xPh;
params.yPh = yPh;

optHypPh = cell(size(Y,2));
parfor nNeuron = 1:size(Y,2)
    optHypPh{nNeuron} = minimize(hypPh, @gp, -100, ...
        @infGaussLik, [], covFuncPh, likFuncPh, xPh, yPh(:,nNeuron));    
end

optHypPh = cat(1,optHypPh{:});

%%
oH = optHyp;
oPh = optHypPh;
gPar = params;
tR = rgTuningCurvesGP(allResp, oH, oPh, gPar);
save(outputDataFilePath,'oH','oPh','gPar','tR'),
end


function tuneResults = rgTuningCurvesGP(allResp, oH, oPh, gPar)

tuneResults = struct;
yTrain = nan(size(gPar.Y));
trainCorr = nan(size(gPar.Y,2),1);
ctTrain = nan(size(gPar.yPh));
ctTrainCorr = nan(size(gPar.yPh,2),1);
parfor nNeur = 1:length(oH)
    if ~isempty(oH(nNeur).cov)
        yTrain(:,nNeur) = gp(oH(nNeur), gPar.infFunc, ...
            gPar.meanFunc, gPar.covFunc, gPar.likFunc, gPar.X, gPar.Y(:,nNeur), gPar.X);
        trainCorr(nNeur) = corr(yTrain(:,nNeur),gPar.Y(:,nNeur));

        ctTrain(:,nNeur) = gp(oPh(nNeur), @infGaussLik, [], ...
            gPar.covFuncPh, gPar.likFuncPh, gPar.xPh, gPar.yPh(:,nNeur), gPar.xPh);
        ctTrainCorr(nNeur) = corr(ctTrain(:,nNeur),gPar.yPh(:,nNeur));
    end
end

tuneResults.yTrain = yTrain;
tuneResults.trainCorr = trainCorr;
tuneResults.ctTrain = ctTrain;
tuneResults.ctTrainCorr = ctTrainCorr;

%% Expand around max response and extract 1-d tunings
expDir = -30:5:30;
expSF = -.5:.1:.5;
expTF = -.5:.1:.5;
expSpd = -.5:.1:.5;

tunDir = linspace(1,360,360);
tunSF = linspace(floor(min(allResp.SF)),ceil(max(allResp.SF)),360);
tunTF = linspace(floor(min(allResp.TF)),ceil(max(allResp.TF)),360);
tunSpd = linspace(0,max(sqrt(allResp.spd)),360);

[~, maxPred] = max(yTrain);
neurDir = nan(length(tunDir),length(oH));
neurSF = nan(length(tunSF),length(oH));
neurTF = nan(length(tunTF),length(oH));
neurSpd = nan(length(tunSpd),length(oH));
varDir = nan(length(tunDir),length(oH));
varSF = nan(length(tunSF),length(oH));
varTF = nan(length(tunTF),length(oH));
varSpd = nan(length(tunSpd),length(oH));
parfor nNeur = 1:length(oH)
    if ~isempty(oH(nNeur).cov)
    indMax = maxPred(nNeur);
    maxDir = allResp.Dir(indMax);
    maxSF = allResp.SF(indMax);
    maxTF = allResp.TF(indMax);
    maxSpd = sqrt(allResp.spd(indMax));
    
    tDir = maxDir + expDir*exp(oH(nNeur).cov(1));
    tSF = maxSF + expSF*exp(oH(nNeur).cov(2));
    tTF = maxTF + expTF*exp(oH(nNeur).cov(3));
    tSpd = maxSpd + expSpd*exp(oH(nNeur).cov(4));
    
    tDir = mod(tDir, 360);
    tSF = max(min(tSF, max(tunSF)),min(tunSF));
    tTF = max(min(tTF, max(tunTF)),min(tunTF));
    tSpd = max(min(tSpd, max(tunSpd)),min(tunSpd));

    [expPred, ~, pDir, pSF, pTF, pSpd] = ...
        gridPredGP(tDir, tSF, tTF, tSpd, oH(nNeur), gPar, gPar.Y(:,nNeur));
    [~, maxExpPred] = max(expPred);
    optDir = pDir(maxExpPred);
    optSF = pSF(maxExpPred);
    optTF = pTF(maxExpPred);
    optSpd = pSpd(maxExpPred);
    
    if abs(optDir-maxDir) == max(abs(expDir))
        warning('Direction Expansion Reached Limit'),
%         keyboard,
    elseif abs(optSF-maxSF) == max(abs(expSF))
        warning('SF Expansion Reached Limit'),
%         keyboard,
    elseif abs(optTF-maxTF) == max(abs(expTF))
        warning('TF Expansion Reached Limit'),
%         keyboard,
    elseif abs(optSpd-maxSpd) == max(abs(expSpd))
        warning('Speed Expansion Reached Limit'),
%         keyboard,
    end
    
    [neurDir(:,nNeur), varDir(:,nNeur)] = gridPredGP(tunDir, optSF, optTF, optSpd, ...
        oH(nNeur), gPar, gPar.Y(:,nNeur));
    [neurSF(:,nNeur),varSF(:,nNeur)] = gridPredGP(optDir, tunSF, optTF, optSpd, ...
        oH(nNeur), gPar, gPar.Y(:,nNeur));
    [neurTF(:,nNeur), varTF(:,nNeur)] = gridPredGP(optDir, optSF, tunTF, optSpd, ...
        oH(nNeur), gPar, gPar.Y(:,nNeur));
    [neurSpd(:,nNeur), varSpd(:,nNeur)] = gridPredGP(optDir, optSF, optTF, tunSpd, ...
        oH(nNeur), gPar, gPar.Y(:,nNeur));
    end
end

tuneResults.tunDir = tunDir;
tuneResults.tunSF = tunSF;
tuneResults.tunTF = tunTF;
tuneResults.tunSpd = tunSpd;
tuneResults.neurTuning = cat(3, neurDir, neurSF, neurTF, neurSpd, ...
    varDir, varSF, varTF, varSpd);

end

function [gridPred, gridVar, pDir, pSF, pTF, pSpd] = gridPredGP(tDir, tSF, tTF, tSpd, oH, gPar, Y)

[pDir, pSF, pTF, pSpd] = ndgrid(tDir, tSF, tTF, tSpd);
pCos = cosd(pDir);
pSin = sind(pDir);
pAll = cat(2,pCos(:), pSin(:), pSF(:), pTF(:), pSpd(:));
[~,~,gridPred, gridVar] = gp(oH, gPar.infFunc, ...
    gPar.meanFunc, gPar.covFunc, gPar.likFunc, gPar.X, Y, pAll);

end