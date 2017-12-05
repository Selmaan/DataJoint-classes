%%
% validBin = v1inf.InfluenceBin & 'bin_id=12';
% validBin = v1inf.Influence & 'exp_date<"2017-07-27"';
validBin = v1inf.Influence & 'exp_date="2017-11-06"';
validStim = v1inf.SelfStim & 'self_stim>5';
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
[iD,iC,rM,rV,sN,sP] = fetchn(((validBin & validFilt) & validStim) - v1inf.Target,...
    'inf_dist','inf_naivecorr','inf_regmu','inf_regvar','inf_shuf_n','inf_shuf_p');
iC(isnan(iC)) = 0;

%%
tmpInd = iD>=25;

tmpDist = iD(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
tmpInf = sN(tmpInd);

dBin = nan(650,1);
for i=1:650
    dBin(i) = nanmean(tmpInf(tmpDist>i & tmpDist<i+50));
end

hold on,
plot(26:675,dBin,'linewidth',2)
axis tight,
xlabel('Neuron-Target Distance (um)')
%%
tmpInd = iD>30;
tmpCorr = iC(tmpInd);
tmpInf = sN(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
% cVals = [-inf,linspace(min(tmpCorr)-.01,max(tmpCorr)+.01,13),inf];
% cVals = linspace(-.4,.7,10);
cVals = [-1, linspace(-.15,.3,7), 1];

cBin = nan(length(cVals)-1,1);
xBin = nan(length(cVals)-1,1);
numInBin = nan(length(cVals)-1,1);
for i=1:length(cVals)-1
    validInds=tmpCorr>=cVals(i) & tmpCorr<cVals(i+1);
    cBin(i) = nanmean(tmpInf(validInds));
    xBin(i) = nanmean(tmpCorr(validInds));
    numInBin(i) = sum((tmpCorr>=cVals(i) & tmpCorr<cVals(i+1)));
end

hold on,plot(xBin,cBin,'linewidth',2)
hold on,plot(xBin,cBin,'*','markersize',10),
% figure,semilogy(xBin,numInBin,'linewidth',2)
% hold on,plot(xBin,numInBin,'*','markersize',10),

%% Distance Only GP Model

hyp = struct();
likFunc = {@likGauss};
hyp.lik = 0;
% covFunc = {@covSEiso};
% hyp.cov = [3;-1];
covFunc = {@covRQiso};
hyp.cov = [3;-1;0];
meanFunc = {@meanConst};
hyp.mean = 0;
infFunc = @infGaussLik;

xu = linspace(25,700,30)'; covFunc = {'apxSparse',covFunc,xu};      % inducing points
infApx  = @(varargin) infFunc(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
hyp.xu = xu;

tmpInd = iD>=30;
tmpDist = iD(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
tmpInf = sN(tmpInd);
optHyp = minimize(hyp,@gp,-60,infApx,meanFunc,covFunc,likFunc,tmpDist,tmpInf);

xDisp = (25:650)';
[~,~,yMu,yVar] = gp(optHyp,infApx,meanFunc,covFunc,likFunc,tmpDist,tmpInf,xDisp);
figure,plot(xDisp,yMu,'k')
hold on,plot(xDisp,yMu+sqrt(yVar),'r')
hold on,plot(xDisp,yMu-sqrt(yVar),'r')

%% Correlation Only GP Model
tmpInd = iD>=30;
tmpCorr = iC(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
% tmpCorr(isnan(tmpInf))=[]; tmpInf(isnan(tmpInf))=[];
tmpInf = sN(tmpInd);

hyp = struct();
likFunc = {@likGauss};
hyp.lik = 0;
% covFunc = {@covSEiso};
% hyp.cov = [-2;-1];
covFunc = {@covRQiso};
hyp.cov = [-2;-1;0];
meanFunc = {@meanConst};
hyp.mean = 0;
infFunc = @infGaussLik;

xu = linspace(min(tmpCorr),max(tmpCorr),20)'; covFunc = {'apxSparse',covFunc,xu};      % inducing points
infApx  = @(varargin) infFunc(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
hyp.xu = xu;

optHyp = minimize(hyp,@gp,-60,infApx,meanFunc,covFunc,likFunc,tmpCorr,tmpInf);

xDisp = linspace(min(tmpCorr),max(tmpCorr),100)';
[~,~,yMu,yVar] = gp(optHyp,infApx,meanFunc,covFunc,likFunc,tmpCorr,tmpInf,xDisp);
figure,plot(xDisp,yMu,'k')
hold on,plot(xDisp,yMu+sqrt(yVar),'r')
hold on,plot(xDisp,yMu-sqrt(yVar),'r')

%% 2-d regression, additive covariance
tmpInd = iD>=30;
tmpDist = iD(tmpInd);
tmpCorr = iC(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
% tmpCorr(isnan(tmpInf))=[]; tmpDist(isnan(tmpInf)) = []; tmpInf(isnan(tmpInf))=[];
tmpInf = sN(tmpInd);

hyp = struct();
likFunc = {@likGauss};
hyp.lik = 0;
meanFunc = {@meanConst};
hyp.mean = 0;
infFunc = @infGaussLik;
% covFun1 = {@covMask, {[1 0], {@covRQ,'iso',[]}}}';
% covFun2 = {@covMask, {[0 1], {@covRQ,'iso',[]}}}';
% covFunc = {@covScale {@covSum, {covFun1, covFun2}}};
% hyp.cov = [3;0;-2;0;0];
covFunc = {@covADD,{[1,2],'covRQiso'}};
hyp.cov = [3;0;0;-2;0;0;0;0];
% covFunc = {@covADD,{[1,2],'covSEiso'}};
% hyp.cov = [3;0;-2;0;0;0];

xu_d = linspace(25,700,10)'; 
xu_c = linspace(min(tmpCorr),max(tmpCorr),10)';
[xu_d,xu_c] = ndgrid(xu_d,xu_c);
xu = [xu_d(:), xu_c(:)];

covFunc = {'apxSparse',covFunc,xu};      % inducing points
infApx  = @(varargin) infFunc(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
hyp.xu = xu;

optHyp = minimize(hyp,@gp,-60,...
    infApx,meanFunc,covFunc,likFunc,[tmpDist,tmpCorr],tmpInf);

[disp_d, disp_c] = ndgrid(linspace(25,650,100)',...
    linspace(min(tmpCorr),max(tmpCorr),99)'); 
xDisp = [disp_d(:), disp_c(:)];
[~,~,yMu,yVar] = gp(optHyp,infApx,meanFunc,covFunc,likFunc,[tmpDist,tmpCorr],tmpInf,xDisp);
yMu = reshape(yMu,100,99); yVar = reshape(yVar,100,99);
figure,imagesc(disp_c(1,:),disp_d(:,1),yMu)
colorbar
xlabel('Neuron-Target Trace Correlation')
ylabel('Neuron Target Distance (um)')

figure,plot(disp_c(1,:),yMu([4,18,50,75,95],:)','linewidth',2)
axis tight
xlabel('Neuron-Target Correlation')
ylabel('Mean Influence')

figure,plot(disp_d(:,1),yMu(:,[5,25,50,75,95]),'linewidth',2)
legend('-.3','-.1','.06','.32','.58')
axis tight
xlabel('Neuron-Target Distance (um)')
ylabel('Mean Influence')
title('Additive Covariance'),

%% Phase-Cluster Analysis

clear
validBin = v1inf.Influence & 'exp_date<"2017-07-27"';
% validBin = v1inf.Influence;
validStim = proj(v1inf.Target & (v1inf.Influence & 'inf_shuf_p>5'));
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
validInf = ((validBin & validFilt) & validStim) - v1inf.Target;
targetTuning = proj(v1inf.Target * v1inf.TuningProps,...
    'neur_id->tmp','test_corr->targ_test_corr','train_corr->targ_train_corr','phase_resp->targ_phase_resp',...
    'pred_resp->targ_pred_resp','resid_resp->targ_resid_resp');
tunedInfluence = validInf * targetTuning * v1inf.TuningProps;

[iD,iC,rM,rV,sN,sP,nTC,tTC,nTRC,tTRC,nPR,tPR,nPred,tPred,nRes,tRes] = fetchn(tunedInfluence,...
    'inf_dist','inf_naivecorr','inf_regmu','inf_regvar','inf_shuf_n','inf_shuf_p',...
    'test_corr','targ_test_corr','train_corr','targ_train_corr','phase_resp','targ_phase_resp',...
    'pred_resp','targ_pred_resp','resid_resp','targ_resid_resp');

sigCorr=nan(length(nPR),1); noiCorr=nan(length(nPR),1); sumCorr=nan(length(nPR),1);
parfor i=1:length(nPR)
    if isempty(nPR{i})
        nPR{i} = nan(120,1);
    end
    if isempty(tPR{i})
        tPR{i} = nan(120,1);
    end
    if ~isempty(nPred{i}) && ~isempty(tPred{i})
        sigCorr(i) = corr(nPred{i},tPred{i});
        noiCorr(i) = corr(nRes{i},tRes{i});
        sumCorr(i) = corr(nPred{i}+nRes{i},tPred{i}+tRes{i});
    end
end
nPR = cat(2,nPR{:});
tPR = cat(2,tPR{:});

load('F:\phaseClusters.mat'),
figure,plot(phaseClusters,'linewidth',1),legend(phaseLabels),
neurDist = nan(length(nTC),size(phaseClusters,2));
targDist = nan(length(tTC),size(phaseClusters,2));
for i=1:size(phaseClusters,2)
    neurDist(:,i) = sqrt(sum((nPR-phaseClusters(:,i)).^2));
    targDist(:,i) = sqrt(sum((tPR-phaseClusters(:,i)).^2));
end
[~,neurIDX] = min(neurDist,[],2);
neurIDX(isnan(sum(neurDist,2))) = nan;
[~,targIDX] = min(targDist,[],2);
targIDX(isnan(sum(targDist,2))) = nan;

tmpInd = iD>=30 & ~isnan(neurIDX) & ~isnan(targIDX);
[p,tbl,stats]=anovan(sN(tmpInd),{neurIDX(tmpInd),targIDX(tmpInd)},'model','interaction','varnames',{'Neuron','Target'});
% [p,tbl,stats]=anovan(rM(tmpInd)./sqrt(rV(tmpInd)),{neurIDX(tmpInd),targIDX(tmpInd)},'model','interaction','varnames',{'Neuron','Target'});
stats.grpnames{1} = phaseLabels(:); stats.grpnames{2} = phaseLabels(:);
figure,[c,m] = multcompare(stats,'Dimension',[1 2]);
figure,multcompare(stats,'Dimension',[1]);
figure,multcompare(stats,'Dimension',[2]);
figure,imagesc(reshape(m(:,1),5,5)),axis square,colorbar,


%% Tuning Curve Analysis
clear
validBin = v1inf.Influence & 'exp_date<"2017-07-27"';
% validBin = v1inf.Influence;
validStim = proj(v1inf.Target & (v1inf.Influence & 'inf_shuf_p>5'));
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
validInf = ((validBin & validFilt) & validStim) - v1inf.Target;

targetTuning = proj(v1inf.Target * v1inf.TuningProps,...
    'neur_id->tmp','test_corr->targ_test_corr','train_corr->targ_train_corr'...
    ,'dir_mean->targ_mdir','dir_std->targ_sdir'...
    ,'sf_mean->targ_msf','sf_std->targ_ssf'...
    ,'tf_mean->targ_mtf','tf_std->targ_stf'...
    ,'spd_mean->targ_mspd','spd_std->targ_sspd'...
    );
tunedInfluence = validInf * targetTuning * v1inf.TuningProps;

[iD,iC,rM,rV,sN,sP,nTC,nTRC,tTC,tTRC,...
    tMD,tSD,nMD,nSD,tMS,tSS,nMS,nSS,tMT,tST,nMT,nST,tMR,tSR,nMR,nSR] = fetchn(tunedInfluence,...
    'inf_dist','inf_naivecorr','inf_regmu','inf_regvar','inf_shuf_n','inf_shuf_p',...
    'test_corr','train_corr','targ_test_corr','targ_train_corr',...
    'targ_mdir','targ_sdir','dir_mean','dir_std',...
    'targ_msf','targ_ssf','sf_mean','sf_std',...
    'targ_mtf','targ_stf','sf_mean','sf_std',...
    'targ_mspd','targ_sspd','spd_mean','spd_std'...
    );

nPairs = length(iD);
norm_tdir = nan(360,nPairs);
norm_ndir = nan(360,nPairs);
norm_tsf = nan(360,nPairs);
norm_nsf = nan(360,nPairs);
norm_ttf = nan(360,nPairs);
norm_ntf = nan(360,nPairs);
norm_tspd = nan(360,nPairs);
norm_nspd = nan(360,nPairs);
for i=1:nPairs
    if ~isempty(tMD{i})
        norm_tdir(:,i) = tMD{i}/mean(tSD{i});
        norm_tsf(:,i) = tMS{i}/mean(tSS{i});
        norm_ttf(:,i) = tMT{i}/mean(tST{i});
        norm_tspd(:,i) = tMR{i}/mean(tSR{i});
    end
    if ~isempty(nMD{i})
        norm_ndir(:,i) = nMD{i}/mean(nSD{i});
        norm_nsf(:,i) = nMS{i}/mean(nSS{i});
        norm_ntf(:,i) = nMT{i}/mean(nST{i});
        norm_nspd(:,i) = nMR{i}/mean(nSR{i});
    end
end
