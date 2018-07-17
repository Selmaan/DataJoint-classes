%{
# Average Residual activity with shuffled target labels
-> v1inf.Influence
-----
avg_rde_shuf: blob      # Average residual de
avg_sv_shuf: blob       # std of residual activity on non-stim trials
avg_shuflabel: blob     # shuffled target label used for defining 'control' sites
%}

classdef ShuffledAvgDeResid < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            if strcmp(fetch1(v1inf.ExpType & key,'exp_type'), 'Multi-Contrast')
                warning('Skipping Multi-Contrast Experiment'),
                return
            end
            targ_label = fetchn(v1inf.Target & key,...
               'targ_label','ORDER BY targ_id');
          
            
            theseStimInfo = v1inf.StimGratingsTraces & key;
            [deResp, stimID, stimVisDir] = fetch1(theseStimInfo,...
                'stimt_de_resp','stimt_targ_id','stimt_vis_dir');
            visDir = stimVisDir(:,1);
            deResp = double(deResp);
            
            binRange = 16:26; %indices for binning response (should be 0-333ms after stim)
            trialResps = squeeze(sum(deResp(binRange,:,:)));

            
            infDist = fetchn(v1inf.Influence & key,'inf_dist','ORDER BY targ_id, neur_id');
            infDist = reshape(infDist,size(deResp,3),[]);

            nShuffles = 1e3;
            allDeMat = nan([size(infDist),nShuffles]);
            allStdVal = nan(size(infDist,1),nShuffles);
            allShufLabels = nan(size(infDist,2),nShuffles);
            parfor nShuf = 1:nShuffles
                perm_order = randperm(length(targ_label));
                shuf_label = targ_label(perm_order);
                [deMat, stdVal] = ...
                    calcAvgResid(trialResps, infDist, stimID, visDir,shuf_label);
                allDeMat(:,:,nShuf) = deMat;
                allStdVal(:,nShuf) = stdVal;
                allShufLabels(:,nShuf) = shuf_label;
            end

            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = nNeuron;
                    keys(nNeuron,nTarg).targ_id = nTarg;
                    keys(nNeuron,nTarg).avg_rde_shuf = squeeze(allDeMat(nNeuron, nTarg,:));
                    keys(nNeuron,nTarg).avg_sv_shuf = squeeze(allStdVal(nNeuron,:));
                    keys(nNeuron,nTarg).avg_shuflabel = squeeze(allShufLabels(nTarg,:));
                end
            end

            insert(self,keys),
        end
    end
end

function [respMat, stdVal] = calcAvgResid(trialResps, infDist, stimID, stimVisDir, targ_label)

distThresh = 30;
[nRespCells,nTargets] = size(infDist);

%% Calculate Residuals
allDirs = unique(stimVisDir);
visAvg = nan(length(allDirs),nRespCells); 
visResid = nan(size(trialResps));
controlStim = find(targ_label <= 0.01);
if isempty(controlStim)
%     warning('No Control Stim in this experiment?'),
    controlStim = find(targ_label<1/2);
end

for s=1:length(allDirs)
    theseTrials = (stimVisDir==allDirs(s));
    for n=1:nRespCells
        validStim = intersect(controlStim, find(infDist(n,:) >=distThresh));
        validInd = ismember(stimID, validStim) & theseTrials;
        visAvg(s,n) = mean(trialResps(validInd,n));
    end
    visResid(theseTrials,:) = trialResps(theseTrials,:) - visAvg(s,:);
end

for n=1:nRespCells
    validStim = ismember(stimID, find(infDist(n,:) >=distThresh));
    stdVal(n) = std(visResid(validStim,n));
end

%% Calculate stim-triggered average traces
respMat = nan(nRespCells,nTargets);
for nStimCell = 1:nTargets
    theseTri = stimID == nStimCell;
    respMat(:,nStimCell) = mean(visResid(theseTri,:));
end

end