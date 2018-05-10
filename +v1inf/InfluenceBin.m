%{
# Influence Measurement Directed Pair
-> v1inf.StimGratingsBinData
-> v1inf.Neuron
-> v1inf.Target
-----
inf_dist: double            # Pairwise distance (in um)
inf_naivecorr=NULL: double  # Neuron-Target Trace Correlation
inf_regmu=NULL: double      # Regression-based estimate of mean effect
inf_regvar=NULL: double     # Regression-based variance of mean estimate
inf_shuf_p=Null: double      # Parametric shuffle estimate of significance
inf_shuf_n=Null: double     # Nonparametric shuffle estimate of significance
%}

classdef InfluenceBin < dj.Computed
    properties
        popRel = v1inf.StimGratingsBinData;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            theseNeurons = v1inf.Neuron & key;
            theseTargets = v1inf.Target & key;
            theseStimInfo = v1inf.StimGratingsBinData & key;
            theseSyncInfo = v1inf.ExpSync & key;
            
            [targ_id, targ_xc, targ_yc, targ_label, targ_neur_id] = fetchn(...
                theseTargets,'targ_id','targ_xc','targ_yc','targ_label','neur_id','ORDER BY targ_id');
            [neur_xc, neur_yc, neur_deconv] = fetchn(...
                theseNeurons,'neur_xc','neur_yc','neur_deconv','ORDER BY neur_id');
            stimData = fetch(theseStimInfo,'*');
            syncData = fetch(theseSyncInfo,'stim_blocks','frame_times');
            
            infDist = sqrt((neur_xc-targ_xc').^2 + (neur_yc-targ_yc').^2);
            naiveCorrs = computeNaiveCorr(neur_deconv, targ_neur_id, syncData);
            if sum(targ_label==0)>0
                [regMu,regVar] = computeInfluenceRegression(stimData,targ_label);
            else
                regMu = nan(size(infDist));
                regVar = nan(size(infDist));
                warning('Experiment on %s contains no control stimulation',key.exp_date),
            end
            [nRespMat, pRespMat] = calcStimShuffle(stimData, infDist, targ_label);
            
            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = nNeuron;
                    keys(nNeuron,nTarg).targ_id = nTarg;
                    keys(nNeuron,nTarg).inf_dist = infDist(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_naivecorr = naiveCorrs(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_regmu = regMu(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_regvar = regVar(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_shuf_p = pRespMat(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_shuf_n = nRespMat(nNeuron,nTarg);
                end
            end
            
            insert(self,keys),
            
        end
    end
end

function naiveCorrs = computeNaiveCorr(neur_deconv, targ_neur_id, syncData)


neur_deconv = cell2mat(neur_deconv');
nNeurons = size(neur_deconv,2);
nTargets = length(targ_neur_id);
validTarg = find(~isnan(targ_neur_id));

validBlocks = find(~syncData.stim_blocks);
naiveCorrs = nan(nNeurons,nTargets,length(validBlocks));
for nBlock = validBlocks
    blockOffsetFrame = length(cat(1,syncData.frame_times{1:nBlock-1}));
    blockInd = blockOffsetFrame + (1:length(syncData.frame_times{nBlock}));
    bin_deconv = squeeze(mean(reshape(neur_deconv(blockInd,:),10,[],nNeurons),1));
    naiveCorrs(:,validTarg,validBlocks==nBlock) = corr(bin_deconv,bin_deconv(:,targ_neur_id(validTarg)));
end
naiveCorrs = nanmean(naiveCorrs,3);
end

function [regMu,regVar] = ...
    computeInfluenceRegression(stimData,targ_label)

% Make all 0 contrast trials have a grating direction of '0'
% Note this works only for multi-contrast with 2 levels where one is zero
if size(stimData.stim_vis_dir,2)==2
    zeroContrasts = stimData.stim_vis_dir(:,2)==0;
    if ~sum(zeroContrasts)
        error('Data Formatting error, see above in code?'),
    end
    stimVisDir = stimData.stim_vis_dir(:,1);
    stimVisDir(zeroContrasts) = 0;
elseif size(stimData.stim_vis_dir,2)==1
    stimVisDir = stimData.stim_vis_dir(:,1);
else
    error('Unknown Data Format'),
end

validStim = find(targ_label>0);
controlStim = find(targ_label == 0);
if isempty(controlStim)
    error('No Control Stimulation Targets Found'),
end
nStim = length(validStim);

allDir = unique(stimVisDir);
X = zeros(size(stimData.stim_de_resp,1),nStim+length(allDir));
for iStim = 1:nStim
    X(:,iStim) = (stimData.stim_targ_id == validStim(iStim));
end

for iDir = 1:length(allDir)
    X(:,nStim+iDir) = (stimVisDir == allDir(iDir));
end

logSpd = log(stimData.stim_mv_spd);
mvBins = 0:25:100;
for mvBin = 1:length(mvBins)-1
    binMin = prctile(logSpd,mvBins(mvBin));
    binMax = prctile(logSpd,mvBins(mvBin+1));
    X(:,nStim+length(allDir)+mvBin) = ...
        (logSpd <= binMax) & (logSpd > binMin);
end
X(:,nStim+length(allDir)+1) = []; % Use lowest movement bin as 'intercept'

validTrials = ~isnan(stimData.stim_targ_id);
X = X(validTrials,:);
yRaw = stimData.stim_de_resp(validTrials,:);

fitMu_raw = nan(size(X,2),size(yRaw,2));
fitVar_raw = nan(size(X,2),size(yRaw,2));
blmObj = diffuseblm(size(X,2),'Intercept',false);
for n=1:size(yRaw,2)
    [~,tmpMu,tmpCov] = estimate(blmObj,X,yRaw(:,n),'Display',false);
    tmpVar = diag(tmpCov);
    fitMu_raw(:,n) = tmpMu;
    fitVar_raw(:,n) = tmpVar;
end

%% Reformat Fitted parameters
regMu = nan(size(yRaw,2),length(targ_label));
regVar = nan(size(yRaw,2),length(targ_label));

regMu(:,validStim) = fitMu_raw(1:nStim,:)';
regVar(:,validStim) = fitVar_raw(1:nStim,:)';
end


function [nRespMat, pRespMat] = ...
    calcStimShuffle(stimData, infDist, targ_label)

distThresh = 30;
nShuffles = 1e5;
[nRespCells,nTargets] = size(infDist);

% Make all 0 contrast trials have a grating direction of '0'
% Note this works only for multi-contrast with 2 levels where one is zero
if size(stimData.stim_vis_dir,2)==2
    zeroContrasts = stimData.stim_vis_dir(:,2)==0;
    if ~sum(zeroContrasts)
        error('Data Formatting error, see above in code?'),
    end
    stimVisDir = stimData.stim_vis_dir(:,1);
    stimVisDir(zeroContrasts) = 0;
elseif size(stimData.stim_vis_dir,2)==1
    stimVisDir = stimData.stim_vis_dir(:,1);
else
    error('Unknown Data Format'),
end


%% Calculate Residual Signal
allDirs = unique(stimVisDir);
visAvg = []; visResid = [];
for s=1:length(allDirs)
    theseTrials = (stimVisDir==allDirs(s));
    for n=1:nRespCells
        validStim = find(infDist(n,:) >= distThresh);
        validInd = theseTrials & ismember(stimData.stim_targ_id,validStim);
        visAvg(s,n) = mean(stimData.stim_de_resp(validInd,n));
    end
    visResid(theseTrials,:) = stimData.stim_de_resp(theseTrials,:) - visAvg(s,:);
end

%% Shuffle Bootstrap

nReps = sum(~isnan(stimData.stim_targ_id))/nTargets;
fprintf('Detected Repetitions of Target Stim: %d \n',nReps),
nTrials = size(visResid,1);

validTrialsCell = cell(0);
for n = 1:nRespCells
    validStim = find(infDist(n,:) >= distThresh);
    validTrialsCell{n} = find(ismember(stimData.stim_targ_id,validStim));
end

shufMat = nan(nShuffles,nRespCells);
parfor n=1:nRespCells
    validTrials = validTrialsCell{n};
    nValid = length(validTrials);
    for nShuffle = 1:nShuffles
        subTri = validTrials(randperm(nValid,nReps));
        shufMat(nShuffle,n) = mean(visResid(subTri,n));
    end
end

shufSTD = mean(abs(shufMat))' * 1.253;

respMat = nan(nRespCells,nTargets);
for nStimCell = 1:nTargets
    theseTri = (stimData.stim_targ_id==nStimCell);
    respMat(:,nStimCell) = mean(visResid(theseTri,:))';
end

controlStim = find(targ_label == 0);
if ~isempty(controlStim)
    for n = 1:nRespCells
        validStim = find(infDist(n,:) >= distThresh);
        validControl = intersect(controlStim,validStim);
        if isempty(validControl)
            controlOffset = 0;
        else
            controlOffset = mean(respMat(n,validControl),2);
        end
        respMat(n,:) = respMat(n,:) - controlOffset;
    end
else
    warning('No Control Targets: Shuffle Offset not computed'),
end

pRespMat = respMat./shufSTD;

for n=1:size(respMat,1)
    for t=1:size(respMat,2)
        gtMat(n,t) = sum(shufMat(:,n)>respMat(n,t));
        ltMat(n,t) = sum(shufMat(:,n)<respMat(n,t));
    end
end

nRespMat = log10((nShuffles-gtMat+1)./(nShuffles-ltMat+1));

end