%{
# Stim-triggered average traces and related info
-> v1inf.Influence
-----
star_df: blob        # Average dF/F
star_de: blob        # Average de
star_ay: blob        # average ay
stdr_df: double      # std of post-synaptic neurons dF/F (for normalizing)
stdr_de: double      # std of post-synaptic neurons de (for normalizing)
stdr_ay: double      # std of post-synaptic neurons AY (for normalizing)
%}

classdef InfTraceResid < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            if strcmp(fetch1(v1inf.ExpType & key,'exp_type'), 'Multi-Contrast')
                warning('Skipping Multi-Contrast Experiment'),
                return
            end
            theseStimInfo = v1inf.StimGratingsTraces & key;
            [dfResp, deResp, ayResp, stimID, stimVisDir] = fetch1(theseStimInfo,...
                'stimt_df_resp','stimt_de_resp','stimt_ay_resp','stimt_targ_id','stimt_vis_dir');
            dfResp = double(dfResp); deResp = double(deResp); ayResp = double(ayResp);
            
            infDist = fetchn(v1inf.Influence & key,'inf_dist','ORDER BY targ_id, neur_id');
            infDist = reshape(infDist,size(dfResp,3),[]);

            [dfMat, dfStd] = calcAvgResid(dfResp, infDist, stimID, stimVisDir(:,1));
            [deMat, deStd] = calcAvgResid(deResp, infDist, stimID, stimVisDir(:,1));
            [ayMat, ayStd] = calcAvgResid(ayResp, infDist, stimID, stimVisDir(:,1));
            
            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = nNeuron;
                    keys(nNeuron,nTarg).targ_id = nTarg;
                    keys(nNeuron,nTarg).stdr_df = dfStd(nNeuron,nTarg);
                    keys(nNeuron,nTarg).stdr_de = deStd(nNeuron,nTarg);
                    keys(nNeuron,nTarg).stdr_ay = ayStd(nNeuron,nTarg);
                    keys(nNeuron,nTarg).star_df = dfMat(:,nNeuron,nTarg);
                    keys(nNeuron,nTarg).star_de = deMat(:,nNeuron,nTarg);
                    keys(nNeuron,nTarg).star_ay = ayMat(:,nNeuron,nTarg);
                end
            end
            
            insert(self,keys),
            
        end
    end
end

function [respMat, stdMat] = calcAvgResid(trialResps, infDist, stimID, stimVisDir)

distThresh = 30;
stdRange = 23:45; %indices for normalization, should be 7-29 frames from stim onset (exclude artifact frames)
[nRespCells,nTargets] = size(infDist);

%% Calculate Residuals and Normalization Factor
allDirs = unique(stimVisDir);
visAvg = nan(size(trialResps,1),length(allDirs),nRespCells); 
visResid = nan(size(trialResps));

for s=1:length(allDirs)
    theseTrials = (stimVisDir==allDirs(s));
    for n=1:nRespCells
        validStim = find(infDist(n,:) >=distThresh);
        validInd = ismember(stimID, validStim) & theseTrials;
        visAvg(:,s,n) = mean(trialResps(:,validInd,n),2);
    end
    visResid(:,theseTrials,:) = trialResps(:,theseTrials,:) - visAvg(:,s,:);
end

for n=1:nRespCells
    validStim = find(infDist(n,:) >=distThresh);
    validResps = visResid(stdRange,validStim,n);
    stdVal(n) = std(mean(validResps));
end

stdMat = repmat(stdVal(:),1,nTargets);

%% Calculate stim-triggered average traces
respMat = nan(size(visResid,1),nRespCells,nTargets);
for nStimCell = 1:nTargets
    theseTri = stimID == nStimCell;
    staTrace = squeeze(mean(visResid(:,theseTri,:),2));
    respMat(:,:,nStimCell) = staTrace;
end

end