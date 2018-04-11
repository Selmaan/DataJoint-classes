%{
# Influence Measurement Directed Pair
-> v1inf.Neuron
-> v1inf.Target
-----
inf_dist: double            # Pairwise distance (in um)
inf_p_cent=Null: double     # Parametric shuffle estimate of significance
inf_n_cent=Null: double     # Nonparametric shuffle estimate of significance
inf_p_raw=Null: double     # Parametric shuffle estimate of significance, without centering on control targets
inf_n_raw=Null: double     # Nonparametric shuffle estimate of significance, without centering on control targets
%}

classdef InfVals < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            theseNeurons = v1inf.Neuron & key;
            theseTargets = v1inf.Target & key;
            theseStimInfo = v1inf.StimGratingsData & key;
            
            [targ_id, targ_xc, targ_yc, targ_label] = fetchn(...
                theseTargets,'targ_id','targ_xc','targ_yc','targ_label','ORDER BY targ_id');
            [neur_id, neur_xc, neur_yc] = fetchn(...
                theseNeurons,'neur_id','neur_xc','neur_yc','ORDER BY neur_id');
            stimData = fetch(theseStimInfo,'*');
            
            infDist = sqrt((neur_xc-targ_xc').^2 + (neur_yc-targ_yc').^2);
            [nRespMat, pRespMat, nRespMat_centered, pRespMat_centered] = calcStimShuffle(stimData, infDist, targ_label);
            
            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = neur_id(nNeuron);
                    keys(nNeuron,nTarg).targ_id = targ_id(nTarg);
                    keys(nNeuron,nTarg).inf_dist = infDist(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_p_raw = pRespMat(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_n_raw = nRespMat(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_p_cent = pRespMat_centered(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_n_cent = nRespMat_centered(nNeuron,nTarg);
                end
            end
            
            insert(self,keys),
            
        end
    end
end

function [nRespMat, pRespMat, nRespMat_centered, pRespMat_centered] = ...
    calcStimShuffle(stimData, infDist, targ_label)

distThresh = 25;
nShuffles = 1e5;
[nRespCells,nTargets] = size(infDist);


%% Calculate Residual Signal
allDirs = unique(stimData.stim_vis_dir);
visAvg = []; visResid = [];
for s=1:length(allDirs)
    theseTrials = (stimData.stim_vis_dir==allDirs(s));
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

% Center data on controls
respMat_centered = respMat;
controlStim = find(targ_label == 0);
for n = 1:nRespCells
    validStim = find(infDist(n,:) >= distThresh);
    validControl = intersect(controlStim,validStim);
    if isempty(validControl)
        warning('No Control Targets for neuron %d',n),
        controlOffset = 0 - mean(shufMat(:,n));
    else
        controlOffset = mean(respMat(n,validControl),2) - mean(shufMat(:,n));
    end
    respMat_centered(n,:) = respMat(n,:) - controlOffset;
end

%%
pRespMat = respMat./shufSTD;
pRespMat_centered = respMat_centered./shufSTD;

for n=1:size(respMat,1)
    for t=1:size(respMat,2)
        gtMat(n,t) = sum(shufMat(:,n)>respMat(n,t));
        ltMat(n,t) = sum(shufMat(:,n)<respMat(n,t));
        gtMat_cent(n,t) = sum(shufMat(:,n)>respMat_centered(n,t));
        ltMat_cent(n,t) = sum(shufMat(:,n)<respMat_centered(n,t));
    end
end

nRespMat = log10(max(nShuffles-gtMat,1)./max(nShuffles-ltMat,1));
nRespMat_centered = log10(max(nShuffles-gtMat_cent,1)./max(nShuffles-ltMat_cent,1));
end