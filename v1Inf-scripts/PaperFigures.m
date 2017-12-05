
clear
% Restrict to Original Dataset minus monitor-off experiments
validExp = (v1inf.Influence - 'exp_date="2017-05-26"' - 'exp_date="2017-05-19"') & 'exp_date<"2017-07-27"';
% Targets must have significant stim responses
validStim = proj(v1inf.Target & (v1inf.Influence & 'inf_shuf_p>5'));
% Neurons and Targets should not overlap
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
% Calculate Overlapping Set
validInf = (validExp & validStim) & validFilt;
% Get tuning similarities
tunedInfluence = validInf * v1inf.PairCorrs;

% Apply only distance criteria
distCriteriaInfluence = validInf & 'inf_dist>25';
% Apply distance, model prediction accuracy and overfitting criteria
tuneCriteriaInfluence = tunedInfluence & 'inf_dist>25' & 'inf_neur_testcorr>0.4' & 'inf_targ_testcorr>0.4' &...
    'inf_neur_traincorr-inf_neur_testcorr < 0.15' & 'inf_targ_traincorr-inf_targ_testcorr < 0.15';

%% Distance Analysis

[iD,traceCorr,sN,sP] = fetchn(distCriteriaInfluence,'inf_dist','inf_naivecorr','inf_shuf_n','inf_shuf_p');
traceCorr(~isfinite(traceCorr)) = 0;

[yM_dist,xM_dist,yS_dist] = valueSmoothingBootstrap(iD,sN,linspace(10,650,250),30,1e3);
figure,plot(xM_dist,yM_dist,'k','linewidth',2)
hold on
plot(xM_dist,yM_dist+yS_dist,'r--','linewidth',2)
plot(xM_dist,yM_dist-yS_dist,'r--','linewidth',2)
axis tight
xlabel('Distance (um)'),

[yM_trace,xM_trace,yS_trace] = valueSmoothingBootstrap(traceCorr,sN,...
    linspace(min(traceCorr),max(traceCorr),250),.1,1e3);
figure,plot(xM_trace,yM_trace,'k','linewidth',2)
hold on
plot(xM_trace,yM_trace+yS_trace,'r--','linewidth',2)
plot(xM_trace,yM_trace-yS_trace,'r--','linewidth',2)
axis tight
xlabel('Trace Correlation'),

figure,plot(xM_trace,yM_trace,'k','linewidth',2)
hold on
plot(xM_trace,yM_trace+yS_trace,'r--','linewidth',2)
plot(xM_trace,yM_trace-yS_trace,'r--','linewidth',2)
axis tight
xlim(prctile(traceCorr,[2.5 97.5])),
xlabel('Trace Correlation'),

%% Tuning Analysis
% Get Pair-Data
[iD,traceCorr,sN,sP,sigCorr,noiCorr,respCorr,dirCorr,oriCorr,sfCorr,tfCorr,spdCorr,tuneCorr] = fetchn(tuneCriteriaInfluence,...
    'inf_dist','inf_naivecorr','inf_shuf_n','inf_shuf_p',...
    'inf_sigcorr','inf_noicorr','inf_respcorr','inf_dircorr',...
    'inf_oricorr','inf_sfcorr','inf_tfcorr','inf_spdcorr',...
    'inf_tuningcorr');

traceCorr(isnan(traceCorr)) = 0;
sigCorr(isnan(sigCorr)) = 0;
noiCorr(isnan(noiCorr)) = 0;
respCorr(isnan(respCorr)) = 0;
dirCorr(isnan(dirCorr))=0;
oriCorr(isnan(oriCorr))=0;
sfCorr(isnan(sfCorr)) = 0;
tfCorr(isnan(tfCorr)) = 0;
spdCorr(isnan(spdCorr)) = 0;
tuneCorr(isnan(tuneCorr)) = 0;

% Make predictor matrices
% Distance slope and offset predictors
dBins=[[25, 100, 300];[100, 300, inf]];
X_dist = nan(length(sN),size(dBins,2)*2);
for dBin = 1:size(dBins,2)
    tmpInd = iD>=dBins(1,dBin) & iD<dBins(2,dBin);
    X_dist(tmpInd,dBin*2) = zscore(iD(tmpInd));
    X_dist(tmpInd,dBin*2 - 1) = 1;
end

% Signal Correlation slope and distance interaction
X_sig = nan(length(sN),2);
X_sig(:,1) = zscore(sigCorr);
X_sig(:,2) = zscore(zscore(iD) .* X_sig(:,1));

% Noise Correlation slope and distance interaction
X_noi = nan(length(sN),2);
X_noi(:,1) = zscore(noiCorr);
X_noi(:,2) = zscore(zscore(iD) .* X_noi(:,1));

% Sig-Noise Interaction
X_int = zscore(X_sig(:,1) .* X_noi(:,1));

%% Dist/Sig/Noi Smoothing Plots

% Signal Corr Plot
xTmp = [X_dist, X_noi]; xTmp(isnan(xTmp)) = 0;
yTmp = sN - xTmp * (pinv(xTmp)*sN);
[yM,xM,yS] = valueSmoothingBootstrap(sigCorr,yTmp,...
    linspace(prctile(sigCorr,0.1),prctile(sigCorr,99.9),250),.1,1e3);
figure,plot(xM,yM,'k','linewidth',2)
hold on
plot(xM,yM+yS,'r--','linewidth',2)
plot(xM,yM-yS,'r--','linewidth',2)
axis tight,
xlim(prctile(sigCorr,[1 99])),
xlabel('Signal Correlation'),

% Noise Corr Plot
xTmp = [X_dist, X_sig]; xTmp(isnan(xTmp)) = 0;
yTmp = sN - xTmp * (pinv(xTmp)*sN);
[yM,xM,yS] = valueSmoothingBootstrap(noiCorr,yTmp,...
    linspace(prctile(noiCorr,0.1),prctile(noiCorr,99.9),250),.1,1e3);
figure,plot(xM,yM,'k','linewidth',2)
hold on
plot(xM,yM+yS,'r--','linewidth',2)
plot(xM,yM-yS,'r--','linewidth',2)
axis tight,
xlim(prctile(noiCorr,[1 99])),
xlabel('Noise Correlation'),

% Trace Corr Plot
xTmp = [X_dist X_sig, X_noi, X_int]; xTmp(isnan(xTmp)) = 0;
yTmp = sN - xTmp * (pinv(xTmp)*sN);
[yM,xM,yS] = valueSmoothingBootstrap(traceCorr,yTmp,...
    linspace(min(traceCorr),max(traceCorr),250),.15,1e3);
figure,plot(xM,yM,'k','linewidth',2)
hold on
plot(xM,yM+yS,'r--','linewidth',2)
plot(xM,yM-yS,'r--','linewidth',2)
axis tight,
xlabel('Trace Correlation'),

figure,plot(xM,yM,'k','linewidth',2)
hold on
plot(xM,yM+yS,'r--','linewidth',2)
plot(xM,yM-yS,'r--','linewidth',2)
axis tight
xlim(prctile(traceCorr,[2.5 97.5])),
xlabel('Trace Correlation'),

%% Shuffle Analysis
nShuffles = 1e4;

xTmp = [X_dist X_sig, X_noi, X_int]; xTmp(isnan(xTmp)) = 0;
xTmp_pinv = pinv(xTmp);
betaShuf = nan(size(xTmp,2),nShuffles);
parfor nShuffle = 1:nShuffles
    betaShuf(:,nShuffle) = xTmp_pinv * sN(randperm(length(sN)));
end

betaTrue = xTmp_pinv * sN;
sigBeta = min(mean(betaTrue<betaShuf,2),mean(betaTrue>betaShuf,2));

%% Extrema Analysis
nShuffles = 1e4;

% Precalculations
isExtrema = traceCorr < prctile(traceCorr,2.5) | traceCorr>prctile(traceCorr,97.5);
xTmp = [X_dist X_sig(:,1), X_noi]; xTmp(isnan(xTmp)) = 0;
xTmp_midRange_pinv = pinv(xTmp(~isExtrema,:));
sN_midRange = sN(~isExtrema);
xTmp_extrema_pinv = pinv(xTmp(isExtrema,:));
sN_extrema = sN(isExtrema);

figure,ecdf(iD(~isExtrema))
hold on,ecdf(iD(isExtrema))

% Fit midRange shuffles
betaShuf = nan(size(xTmp,2),nShuffles);
parfor nShuffle = 1:nShuffles
    betaShuf(:,nShuffle) = xTmp_midRange_pinv * sN_midRange(randperm(length(sN_midRange)));
end
betaTrue_midRange = xTmp_midRange_pinv * sN_midRange;
sigBeta_midRange = min(mean(betaTrue_midRange<betaShuf,2),mean(betaTrue_midRange>betaShuf,2));

% Fit Extrema shuffles
betaShuf = nan(size(xTmp,2),nShuffles);
parfor nShuffle = 1:nShuffles
    betaShuf(:,nShuffle) = xTmp_extrema_pinv * sN_extrema(randperm(length(sN_extrema)));
end
betaTrue_extrema = xTmp_extrema_pinv * sN_extrema;
sigBeta_extrema = min(mean(betaTrue_extrema<betaShuf,2),mean(betaTrue_extrema>betaShuf,2));

figure,bar([betaTrue_midRange(7:end),betaTrue_extrema(7:end)])
[sigBeta_midRange(7:end),sigBeta_extrema(7:end)],