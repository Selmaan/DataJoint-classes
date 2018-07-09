%{
# Average response-normed residual activity (alternative to 'influence')
-> v1inf.Influence
-----
avg_rde: double    # Average response-normed residual de
%}

classdef AvgDeResidRN < dj.Computed
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
            deResp = double(deResp);
            
            infDist = fetchn(v1inf.Influence & key,'inf_dist','ORDER BY targ_id, neur_id');
            infDist = reshape(infDist,size(deResp,3),[]);

            deMat = calcAvgResid(deResp, infDist, stimID, stimVisDir(:,1),targ_label);

            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = nNeuron;
                    keys(nNeuron,nTarg).targ_id = nTarg;
                    keys(nNeuron,nTarg).avg_rde = deMat(nNeuron, nTarg);
                end
            end
            
            insert(self,keys),
            
        end
    end
end

function respMat = calcAvgResid(trialResps, infDist, stimID, stimVisDir, targ_label)

distThresh = 30;
binRange = 16:26; %indices for binning response (should be 0-333ms after stim)

[nRespCells,nTargets] = size(infDist);
trialResps = squeeze(sum(trialResps(binRange,:,:)));

%% Calculate Residuals
allDirs = unique(stimVisDir);
visAvg = nan(length(allDirs),nRespCells); 
visResid = nan(size(trialResps));
controlStim = find(targ_label <= 0.01);
if isempty(controlStim)
    warning('No Control Stim in this experiment?'),
    controlStim = find(targ_label<1/2);
end

for s=1:length(allDirs)
    theseTrials = (stimVisDir==allDirs(s));
    for n=1:nRespCells
        validStim = intersect(controlStim, find(infDist(n,:) >=distThresh));
        validInd = ismember(stimID, validStim) & theseTrials;
        visAvg(s,n) = mean(trialResps(validInd,n));
    end
    visResid(theseTrials,:) = (trialResps(theseTrials,:) - visAvg(s,:))./visAvg(s,:);
end

%% Calculate stim-triggered average traces
respMat = nan(nRespCells,nTargets);
for nStimCell = 1:nTargets
    theseTri = stimID == nStimCell;
    respMat(:,nStimCell) = mean(visResid(theseTri,:));
end

end