%{
# Influence Measurement Directed Pair
-> v1inf.Neuron
-> v1inf.Target
-----
inf_p_0=Null: double     # Parametric shuffle estimate of significance, 0% contrast trials
inf_n_0=Null: double     # Nonparametric shuffle estimate of significance, 0% contrast trials
inf_p_10=Null: double     # Parametric shuffle estimate of significance, 10% contrast trials
inf_n_10=Null: double     # Nonparametric shuffle estimate of significance, 0% contrast trials
%}

classdef MultiInf < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            % Only process multi-contrast experiments
            thisExptType = fetch1(v1inf.ExpType & key,'exp_type');
            if ~strcmp(thisExptType,'Multi-Contrast')
                fprintf('This Experiment not Multi-Contrast \n'),
                return
            end
            
            theseNeurons = v1inf.Neuron & key;
            theseTargets = v1inf.Target & key;
            theseStimInfo = v1inf.StimGratingsData & key;
            
            [targ_id, targ_xc, targ_yc, targ_label] = fetchn(...
                theseTargets,'targ_id','targ_xc','targ_yc','targ_label','ORDER BY targ_id');
            [neur_id, neur_xc, neur_yc] = fetchn(...
                theseNeurons,'neur_id','neur_xc','neur_yc','ORDER BY neur_id');
            stimData = fetch(theseStimInfo,'*');
            
            infDist = sqrt((neur_xc-targ_xc').^2 + (neur_yc-targ_yc').^2);
            
            allContrasts = unique(stimData.stim_vis_dir(:,2));
            if ~all(allContrasts == [0; .1])
                error('Contrasts do not match expected values'),
            end
            t0 = stimData.stim_vis_dir(:,2)==0;
            t10 = stimData.stim_vis_dir(:,2)==0.1;
            stimData_0 = subsetStimTrials(stimData,t0);
            stimData_0.stim_vis_dir(:,1) = 0; % Set all zero-contrast trials to have same 'direction'
            stimData_10 = subsetStimTrials(stimData,t10);
            [nRespMat_0, pRespMat_0] = calcStimShuffle(stimData_0, infDist, targ_label);
            [nRespMat_10, pRespMat_10] = calcStimShuffle(stimData_10, infDist, targ_label);

            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = neur_id(nNeuron);
                    keys(nNeuron,nTarg).targ_id = targ_id(nTarg);
                    keys(nNeuron,nTarg).inf_p_0 = pRespMat_0(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_n_0 = nRespMat_0(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_p_10 = pRespMat_10(nNeuron,nTarg);
                    keys(nNeuron,nTarg).inf_n_10 = nRespMat_10(nNeuron,nTarg);
                end
            end
            
            insert(self,keys),
            
        end
    end
end

function stimData_ind = subsetStimTrials(stimData,tInd)

stimData_ind = stimData;
stimData_ind.stim_targ_id = stimData_ind.stim_targ_id(tInd,1);
stimData_ind.stim_vis_dir = stimData_ind.stim_vis_dir(tInd,1);
stimData_ind.stim_mv_spd = stimData_ind.stim_mv_spd(tInd,1);
stimData_ind.stim_pre_spd = stimData_ind.stim_pre_spd(tInd,1);

stimData_ind.stim_de_resp = stimData_ind.stim_de_resp(tInd,:);
stimData_ind.stim_de_pre = stimData_ind.stim_de_pre(tInd,:);

end

function [nRespMat, pRespMat] = ...
    calcStimShuffle(stimData, infDist, targ_label)

distThresh = 25;
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

nReps = round(sum(~isnan(stimData.stim_targ_id))/nTargets);
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
controlStim = find(targ_label == 0);
for n = 1:nRespCells
    validStim = find(infDist(n,:) >= distThresh);
    validControl = intersect(controlStim,validStim);
    if isempty(validControl)
        warning('No Control Targets for neuron %d',n),
        controlOffset = 0;
    else
        controlOffset = mean(respMat(n,validControl),2);
    end
    respMat(n,:) = respMat(n,:) - controlOffset;
end

%%
pRespMat = respMat./shufSTD;

for n=1:size(respMat,1)
    for t=1:size(respMat,2)
        gtMat(n,t) = sum(shufMat(:,n)>respMat(n,t));
        ltMat(n,t) = sum(shufMat(:,n)<respMat(n,t));
    end
end

nRespMat = log10(max(nShuffles-gtMat,1)./max(nShuffles-ltMat,1));
end