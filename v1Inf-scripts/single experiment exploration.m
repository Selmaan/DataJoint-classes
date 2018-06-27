%% Obtain list of Experiments
clear,
allExpts = fetch(v1inf.ExpType & 'exp_type="Multi-Contrast"','ORDER BY exp_date');
% thisInf = v1inf.Influence * v1inf.PairCorrs;
thisInf = v1inf.MultiInf * v1inf.Influence * v1inf.PairCorrs; % Remember to change below code too!

% Targets must have significant stim responses
validStim = proj(v1inf.Target & (v1inf.Influence & 'inf_shuf_p>5'));
% Neurons and Targets should not overlap
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
% Calculate Overlapping Set
validInf = (thisInf & validStim) & validFilt;

tuneCriteriaInfluence = validInf & 'inf_dist>25' & 'inf_neur_testcorr>0.4' & 'inf_targ_testcorr>0.4' &...
    'inf_neur_traincorr-inf_neur_testcorr < 0.15' & 'inf_targ_traincorr-inf_targ_testcorr < 0.15';

%% Load influence and tuning data for indivual experiment
allBeta = nan(length(allExpts),10);
for nX = 1:length(allExpts) 
    % Get Pair-Data
    [iD,sN,sigCorr,noiCorr] = ...
        fetchn(tuneCriteriaInfluence & allExpts(nX),...
        'inf_dist','inf_n_10','inf_sigcorr','inf_noicorr');
%         'inf_dist','inf_shuf_n','inf_sigcorr','inf_noicorr');
    sigCorr(isnan(sigCorr)) = 0;
    noiCorr(isnan(noiCorr)) = 0;
    
    % Make Predictor Matrix
    dBins=[[25, 125, 300];[125, 300, inf]];
    X_dist = nan(length(sN),size(dBins,2)*2);
    for dBin = 1:size(dBins,2)
        tmpInd = iD>=dBins(1,dBin) & iD<dBins(2,dBin);
        X_dist(tmpInd,dBin*2) = zscore(iD(tmpInd));
        X_dist(tmpInd,dBin*2 - 1) = 1;
    end
    allDist = zscore(log(iD));
    X_sig = nan(length(sN),2);
    X_sig(:,1) = zscore(sigCorr);
    X_sig(:,2) = zscore(allDist .* X_sig(:,1));
    X_noi = nan(length(sN),2);
    X_noi(:,1) = zscore(noiCorr);
    X_noi(:,2) = zscore(allDist .* X_noi(:,1));
    
    % Fit Regression Model
    xTmp = [X_dist X_sig, X_noi]; xTmp(isnan(xTmp)) = 0;
    allBeta(nX,:) = xTmp \ sN;
    
end

sigCoefs = allBeta(:,7);
figure, plot(sigCoefs,'.','markersize',30),

%% Analyze collection of model fits over experiments

% distanceEffects = allBeta(:,1:6);
% figure,boxplot(distanceEffects,...
% 'Colors',.15 * [1 1 1],'MedianStyle','line','BoxStyle','filled')
% a = get(get(gca,'children'),'children');
% t = get(a,'tag');
% idx = strcmp(t,'Box');
% boxes = a(idx);
% set(boxes,'linewidth',15);
% set(a(strcmp(t,'Median')),'Color',[.85 .33 .1]),
% set(a(strcmp(t,'Median')),'linewidth',1.5),
% hold on,
% line([0 6]+1/2,[0 0],'color',1/2 * [1 1 1],'linestyle','--'),
% ylabel('Beta Coef'),
% xticklabels({'Offset (25-125um)','Slope (25-125um)',...
%     'Offset (125-300um)','Slope (125-300um)',...
%     'Offset (>300um)','Slope (>300um)',}),
% xtickangle(45),
% title('Distance Predictors'),
% 
% tuningEffects = allBeta(:,7:end);
% figure,boxplot(tuningEffects,...
% 'Colors',.15 * [1 1 1],'MedianStyle','line','BoxStyle','filled')
% a = get(get(gca,'children'),'children');
% t = get(a,'tag');
% idx = strcmp(t,'Box');
% boxes = a(idx);
% set(boxes,'linewidth',15);
% set(a(strcmp(t,'Median')),'Color',[.85 .33 .1]),
% set(a(strcmp(t,'Median')),'linewidth',1.5),
% hold on,
% line([0 6]+1/2,[0 0],'color',1/2 * [1 1 1],'linestyle','--'),
% ylabel('Beta Coef'),
% xticklabels({'Signal Corr','Sig*Dist',...
%     'Noise Corr','Noise*Dist'}),
% xtickangle(45),
% title('Tuning Predictors'),

%%

%% Use a single LME w/ session-level random effects instead
clear,
allExpts = v1inf.ExpType & 'exp_type="Multi-Contrast"';
% thisInf = v1inf.Influence * v1inf.PairCorrs;
thisInf = v1inf.MultiInf * v1inf.Influence * v1inf.PairCorrs; % Remember to change below code too!

% Targets must have significant stim responses
validStim = proj(v1inf.Target & (v1inf.Influence & 'inf_shuf_p>5'));
% Neurons and Targets should not overlap
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
% Calculate Overlapping Set
validInf = (thisInf & validStim) & validFilt;

tuneCriteriaInfluence = validInf & 'inf_dist>25' & 'inf_neur_testcorr>0.4' & 'inf_targ_testcorr>0.4' &...
    'inf_neur_traincorr-inf_neur_testcorr < 0.15' & 'inf_targ_traincorr-inf_targ_testcorr < 0.15';

[expDate,iD,sN,sigCorr,noiCorr] = ...
    fetchn(tuneCriteriaInfluence & allExpts,'exp_date',...
    'inf_dist','inf_n_10','inf_sigcorr','inf_noicorr','ORDER BY exp_date');

% Make session ID index
allSesh = unique(expDate);
seshID = nan(length(expDate),1);
for iSesh = 1:length(allSesh)
    seshInd = strcmp(expDate,allSesh(iSesh));
    seshID(seshInd) = iSesh;
end
seshID = categorical(seshID);

% Make Predictor Matrix
dBins=[[25, 125, 300];[125, 300, inf]];
X_dist = zeros(length(sN),size(dBins,2)*2);
for dBin = 1:size(dBins,2)
    tmpInd = iD>=dBins(1,dBin) & iD<dBins(2,dBin);
    X_dist(tmpInd,dBin*2) = zscore(iD(tmpInd));
    X_dist(tmpInd,dBin*2 - 1) = 1;
end
allDist = zscore(log(iD));
X_sig = nan(length(sN),2);
X_sig(:,1) = zscore(sigCorr);
X_sig(:,2) = zscore(allDist .* X_sig(:,1));
X_noi = nan(length(sN),2);
X_noi(:,1) = zscore(noiCorr);
X_noi(:,2) = zscore(allDist .* X_noi(:,1));

X = table(X_dist(:,1),X_dist(:,2),X_dist(:,3),X_dist(:,4),X_dist(:,5),X_dist(:,6),...
    X_sig(:,1), X_sig(:,2), X_noi(:,1), X_noi(:,2), seshID, sN);
formula = ['sN ~ -1 + Var1 + Var2 + Var3 + Var4 + Var5 + Var6 + Var7 + Var8 + Var9 + Var10'...
    ' + (1|seshID) + (Var7-1|seshID)'];
mdl = fitlme(X, formula,'DummyVarCoding','effects');
[tmpR, tmpNames] = mdl.randomEffects; tmpF = mdl.fixedEffects;
sigCoefs = tmpR(length(allSesh)+1:end) + tmpF(7);

%% RESULTS OF ABOVE ANALYSIS
% Saved in supplemental folder on dropbox

%% Number of days since injection

clear,
load('sigCoefs_multiContrast.mat'),
[doi, allExpts] = fetchn(v1inf.Mouse * v1inf.ExpType & 'exp_type="Multi-Contrast"',...
    'doi','exp_date','ORDER BY exp_date');
elapsedDays = datenum(allExpts)-datenum(doi);
figure,plot(elapsedDays,sigCoefs_lm,'.','markersize',30),
figure,plot(elapsedDays,sigCoefs_lme,'.','markersize',30),
[c,p] = corr([elapsedDays,sigCoefs_lm,sigCoefs_lme],'type','Spearman'),

%% Number of experiments in the same mouse, and FOV depth

clear,
load('sigCoefs_multiContrast.mat'),
[mID,expDate, zD] = fetchn(v1inf.ExpDepth * v1inf.ExpType & 'exp_type="Multi-Contrast"',...
    'mouse_id','exp_date','fov_depth','ORDER BY exp_date');
allMice = unique(mID);
expNum = nan(length(mID),1);
for n = 1:length(allMice)
    thisMouse = allMice(n);
    mInd = mID==thisMouse;
    expNum(mInd) = 1:sum(mInd);
end
figure,plot(expNum,sigCoefs_lm,'.','markersize',30),
xlabel('# of Experiments in Same Mouse'),ylabel('b-sig Ind. Model'),
figure,plot(expNum,sigCoefs_lme,'.','markersize',30),
xlabel('# of Experiments in Same Mouse'),ylabel('b-sig LME'),
[c,p] = corr([expNum,sigCoefs_lm,sigCoefs_lme],'type','Spearman','rows','complete'),

% Example code to restrict by these
% theseExpt = expDate(expNum<=5 & mID~=37);
% key = struct('exp_date',theseExpt);
%% Running speed, or fraction of time running, within session

clear,
load('sigCoefs_multiContrast.mat'),
allExpts = fetch(v1inf.ExpType & 'exp_type="Multi-Contrast"',...
    'exp_date','ORDER BY exp_date');
runStats = nan(length(allExpts),4);
for n=1:length(allExpts)
    tmp = fetch1(v1inf.StimGratingsData & allExpts(n),'stim_mv_spd');
    runStats(n,:) = [mean(tmp>.15), median(tmp), mean(tmp), mad(tmp)];
end

figure,plot(runStats(:,1),sigCoefs_lm,'.','markersize',30),
xlabel('Fraction of Trials Running'),ylabel('b-sig Ind. Model'),
figure,plot(runStats(:,1),sigCoefs_lme,'.','markersize',30),
xlabel('Fraction of Trials Running'),ylabel('b-sig LME'),
[c,p] = corr([runStats,sigCoefs_lm,sigCoefs_lme],'type','Spearman'),

%% Injection Type

clear,
load('sigCoefs_singleContrast.mat'),
[chr, allExpts] = fetchn(v1inf.Mouse * v1inf.ExpType & 'exp_type="Single-Contrast"',...
    'inj_chr','exp_date','ORDER BY exp_date');
isC1V1 = strcmp(chr,'C1V1');
figure,hold on
plot(sigCoefs_lm(isC1V1),sigCoefs_lme(isC1V1),'.','markersize',30),
plot(sigCoefs_lm(~isC1V1),sigCoefs_lme(~isC1V1),'.','markersize',30),
legend('C1V1','ChrimsonR'),
axis equal,

[c,p] = corr([isC1V1,sigCoefs_lm,sigCoefs_lme],'type','Spearman'),
