%{
# Influence Measurement Directed Pair (Regression-based)
-> v1inf.Neuron
-> v1inf.Target
-----
inf_dist: double            # Pairwise distance (in um)
inf_naivecorr=NULL: double  # Neuron-Target Trace Correlation
inf_rawmu=NULL: double      # Mean influence of stim period response
inf_rawvar=NULL: double     # Variance of stim period response estimate
inf_difmu=Null: double      # Mean influence of stim minus pre-period response
inf_difvar=Null: double     # Variance of stim minus pre-period estimate
%}

classdef Influence < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            theseNeurons = v1inf.Neuron & key;
            theseTargets = v1inf.Target & key;
            theseStimInfo = v1inf.StimGratingsData & key;
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
                [rawMu,rawVar, difMu, difVar] = computeInfluenceRegression(stimData,targ_label);
            else
                rawMu = nan(size(infDist));
                rawVar = nan(size(infDist));
                difMu = nan(size(infDist));
                difVar = nan(size(infDist));
                warning('Experiment on %s contains no control stimulation',key.exp_date),
            end
            
            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = nNeuron;
                    keys(nNeuron,nTarg).targ_id = nTarg;
                    keys(nNeuron,nTarg).inf_dist = infDist(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_naivecorr = naiveCorrs(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_rawmu = rawMu(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_rawvar = rawVar(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_difmu = difMu(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_difvar = difVar(nNeuron,nTarg);
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

function [rawMu,rawVar, difMu, difVar] = ...
    computeInfluenceRegression(stimData,targ_label)

validStim = find(targ_label>0);
controlStim = find(targ_label == 0);
if isempty(controlStim)
    error('No Control Stimulation Targets Found'),
end
nStim = length(validStim);

allDir = unique(stimData.stim_vis_dir);
X = zeros(size(stimData.stim_de_resp,1),nStim+length(allDir));
for iStim = 1:nStim
    X(:,iStim) = (stimData.stim_targ_id == validStim(iStim));
end

for iDir = 1:length(allDir)
    X(:,nStim+iDir) = (stimData.stim_vis_dir == allDir(iDir));
end

mvBins = 0:25:100;
for mvBin = 1:length(mvBins)-1
    binMin = prctile(stimData.stim_mv_spd,mvBins(mvBin));
    binMax = prctile(stimData.stim_mv_spd,mvBins(mvBin+1));
    X(:,nStim+length(allDir)+mvBin) = ...
        (stimData.stim_mv_spd <= binMax) & (stimData.stim_mv_spd > binMin);
end
X(:,nStim+length(allDir)+1) = []; % Use lowest movement bin as 'intercept'

yRaw = stimData.stim_de_resp;
yDif = stimData.stim_de_resp-stimData.stim_de_pre;

fitMu_raw = nan(size(X,2),size(yRaw,2));
fitVar_raw = nan(size(X,2),size(yRaw,2));
fitMu_dif = nan(size(X,2),size(yRaw,2));
fitVar_dif = nan(size(X,2),size(yRaw,2));
blmObj = diffuseblm(size(X,2),'Intercept',false);
for n=1:size(yRaw,2)
    [~,tmpMu,tmpCov] = estimate(blmObj,X,yRaw(:,n),'Display',false);
    tmpVar = diag(tmpCov);
    fitMu_raw(:,n) = tmpMu;
    fitVar_raw(:,n) = tmpVar;
    
    [~,tmpMu,tmpCov] = estimate(blmObj,X,yDif(:,n),'Display',false);
    tmpVar = diag(tmpCov);
    fitMu_dif(:,n) = tmpMu;
    fitVar_dif(:,n) = tmpVar;
end

%% Reformat Fitted parameters
rawMu = nan(size(yRaw,2),length(targ_label));
rawVar = nan(size(yRaw,2),length(targ_label));
difMu = nan(size(yRaw,2),length(targ_label));
difVar = nan(size(yRaw,2),length(targ_label));

rawMu(:,validStim) = fitMu_raw(1:nStim,:)';
rawVar(:,validStim) = fitVar_raw(1:nStim,:)';
difMu(:,validStim) = fitMu_dif(1:nStim,:)';
difVar(:,validStim) = fitVar_dif(1:nStim,:)';

% respBeta = fitMu(nStim+1:end,:);
% respBetaVar = fitVar(nStim+1:end,:);
end