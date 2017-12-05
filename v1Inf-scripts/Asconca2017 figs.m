%% (0) DataJoint Restriction and Table Construction

clear
% Restrict to Original Dataset
validBin = v1inf.Influence & 'exp_date<"2017-07-27"';
% Targets must have significant stim responses
validStim = proj(v1inf.Target & (v1inf.Influence & 'inf_shuf_p>5'));
% Neurons and Targets should not overlap
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
% Calculate Overlapping Set
validInf = ((validBin & validFilt) & validStim) - v1inf.Target;

% Construct Tuning properties table for targets, removing primary key neur_id
targetTuning = proj(v1inf.Target * v1inf.TuningProps,'neur_id->tmp',...
    'test_corr->targ_test_corr','train_corr->targ_train_corr',...
    'pred_resp->targ_pred_resp','resid_resp->targ_resid_resp'...
    ,'dir_mean->targ_mdir','dir_std->targ_sdir'...
    ,'sf_mean->targ_msf','sf_std->targ_ssf'...
    ,'tf_mean->targ_mtf','tf_std->targ_stf'...
    ,'spd_mean->targ_mspd','spd_std->targ_sspd');

% Combine influence with neuron and target tuning
tunedInfluence = validInf * targetTuning * v1inf.TuningProps;
% Apply distance, model prediction accuracy and overfitting criteria
criteriaInfluence = tunedInfluence & 'inf_dist>25' & 'test_corr>0.4' & 'targ_test_corr>0.4' &...
    'train_corr-test_corr < 0.15' & 'targ_train_corr-targ_test_corr < 0.15';
% Fetch from datajoint
[iD,iC,rM,rV,sN,sP,nTC,tTC,nTRC,tTRC,nPred,tPred,nRes,tRes,tMD,tMS,tMT,tMR,nMD,nMS,nMT,nMR] = fetchn(criteriaInfluence,...
    'inf_dist','inf_naivecorr','inf_regmu','inf_regvar','inf_shuf_n','inf_shuf_p',...
    'test_corr','targ_test_corr','train_corr','targ_train_corr',...
    'pred_resp','targ_pred_resp','resid_resp','targ_resid_resp',...
    'targ_mdir','targ_msf','targ_mtf','targ_mspd',...
    'dir_mean','sf_mean','tf_mean','spd_mean');
% Compute signal, noise, and response correlation
sigCorr=nan(length(nTC),1); noiCorr=nan(length(nTC),1); respCorr=nan(length(nTC),1);
dirCorr=nan(length(nTC),1); sfCorr=nan(length(nTC),1); tfCorr=nan(length(nTC),1); spdCorr=nan(length(nTC),1);
oriCorr = nan(length(nTC),1); tuningCorr = nan(length(nTC),1);
parfor i=1:length(nTC)
        sigCorr(i) = corr(nPred{i},tPred{i});
        noiCorr(i) = corr(nRes{i},tRes{i});
        respCorr(i) = corr(nPred{i}+nRes{i},tPred{i}+tRes{i});
        dirCorr(i) = corr(tMD{i},nMD{i});
        oriCorr(i) = corr(tMD{i}(1:180) + tMD{i}(181:360) , nMD{i}(1:180)+nMD{i}(181:360));
        sfCorr(i) = corr(tMS{i},nMS{i});
        tfCorr(i) = corr(tMT{i},nMT{i});
        spdCorr(i) = corr(tMR{i},nMR{i});
        tuningCorr(i) = corr([tMD{i};tMS{i};tMT{i};tMR{i}],[nMD{i};nMS{i};nMT{i};nMR{i}]);
end
dirCorr(isnan(dirCorr))=0;
oriCorr(isnan(oriCorr))=0;
sfCorr(isnan(sfCorr))=0;
tfCorr(isnan(tfCorr))=0;
spdCorr(isnan(spdCorr))=0;

%% (1) GP models

gpDist = infV1_gp1d(iD,sN);
gpTraceCorr = infV1_gp1d(iC,sN);
gpSigCorr = infV1_gp1d(sigCorr,sN);

for tmp = [gpDist,gpSigCorr,gpTraceCorr]
    figure,plot(tmp.xDisp,tmp.yMu,'k','linewidth',2)
    hold on,plot(tmp.xDisp,tmp.yMu+tmp.yStd,'r--','linewidth',2)
    hold on,plot(tmp.xDisp,tmp.yMu-tmp.yStd,'r--','linewidth',2)
    axis tight,
%     xlim([20 600]),
end

gpDistRC = infV1_gp2d_lin([iD,respCorr],sN);
yMu = reshape(gpDistRC.yMu,100,100);
xDisp = reshape(gpDistRC.xDisp,100,100,2);
figure,plot(xDisp(:,[1 100 50],1),yMu(:,[1 100 50]),'linewidth',2)
axis tight, xlim([20 600]),

%% (2) Piece-wise regression

dBins=[[25, 100, 300];[100, 300, inf]];
nShuffles = 1e3;
opt = glmnetSet(struct('alpha',1e-3,'nlambda',250,'lambda_min',1e-5));

rawBeta = nan(3,3);
normBeta = nan(3,3);
sigBeta = nan(3,3);

for dBin = 1:size(dBins,2)
    dBin,
    tmpInd = iD>dBins(1,dBin) & iD<dBins(2,dBin);
    X = [iD(tmpInd),sigCorr(tmpInd),noiCorr(tmpInd)];
    Yn = sN(tmpInd);
    
    bShuf = nan(size(X,2),nShuffles);
    parfor nShuffle = 1:nShuffles
        shufOrder = randperm(size(X,1));
        data = cvglmnet(X(shufOrder,:),Yn,'gaussian',opt,'deviance',20);
        bShuf(:,nShuffle) = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
    end
    data = cvglmnet(X,Yn,'gaussian',opt,'deviance',20);
    bOpt = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
    rawBeta(dBin,:) = bOpt;
    normBeta(dBin,:) = bOpt.*nanstd(X)';
    sigBeta(dBin,:) = mean(bOpt > bShuf,2);
end

figure,bar(normBeta),
xticklabels({'25-100um','100-300um','>300 um'}),
ylim([-.085 .06])
yticks(-.08:.02:.06)
xlabel('Distance Bin')
ylabel('Influence Slope')
legend('Distance','Signal Corr','Noise Corr')

figure,scatter(sigCorr,noiCorr,6,iD,'filled')
axis tight,
xlabel('Signal Correlation')
ylabel('Noise Correlation')

% gpDistCorrS = infV1_gp1d(iD,sigCorr);
% gpDistCorrN = infV1_gp1d(iD,noiCorr);
% figure,hold on,
% plot(gpDistCorrS.xDisp,gpDistCorrS.yMu,'b','linewidth',2),
% plot(gpDistCorrS.xDisp,gpDistCorrS.yMu+gpDistCorrS.yStd,'b--','linewidth',2),
% plot(gpDistCorrS.xDisp,gpDistCorrS.yMu-gpDistCorrS.yStd,'b--','linewidth',2),
% plot(gpDistCorrN.xDisp,gpDistCorrN.yMu,'r','linewidth',2),
% plot(gpDistCorrN.xDisp,gpDistCorrN.yMu+gpDistCorrN.yStd,'r--','linewidth',2),
% plot(gpDistCorrN.xDisp,gpDistCorrN.yMu-gpDistCorrN.yStd,'r--','linewidth',2),


%% (3) Sub-tuning with piece-wise regression

nShuffles = 1e3;
allTuneVar = [sigCorr,tuningCorr,dirCorr,oriCorr,sfCorr,tfCorr,spdCorr];
nTuneVar = size(allTuneVar,2);
opt = glmnetSet(struct('alpha',1e-3,'nlambda',250,'lambda_min',1e-5));

rawBeta = nan(nTuneVar,3);
normBeta = nan(nTuneVar,3);
sigBeta = nan(nTuneVar,3);
for xTune = 1:nTuneVar
    xTune,
    X = [iD,allTuneVar(:,xTune),noiCorr];
    Yn = sN;

    bShuf = nan(size(X,2),nShuffles);
    parfor nShuffle = 1:nShuffles
        shufOrder = randperm(size(X,1));
        data = cvglmnet(X(shufOrder,:),Yn,'gaussian',opt,'deviance',20);
        bShuf(:,nShuffle) = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
    end
    data = cvglmnet(X,Yn,'gaussian',opt,'deviance',20);
    bOpt = data.glmnet_fit.beta(:,data.lambda==data.lambda_min);
    rawBeta(xTune,:) = bOpt;
    normBeta(xTune,:) = bOpt.*nanstd(X)';
    sigBeta(xTune,:) = mean(bOpt > bShuf,2);
end

figure,imagesc(corrcoef([iD,noiCorr,allTuneVar]),[-.05 1]),
figure,bar(normBeta),