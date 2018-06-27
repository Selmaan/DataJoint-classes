%{
# Stim-triggered average traces and related info
-> v1inf.Influence
-----
sta_df: blob        # Average dF/F
sta_de: blob        # Average de
sta_ay: blob        # average ay
std_df: double      # std of post-synaptic neurons dF/F (for normalizing)
std_de: double      # std of post-synaptic neurons de (for normalizing)
std_ay: double      # std of post-synaptic neurons AY (for normalizing)
%}

classdef InfTrace < dj.Computed
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
            [dfResp, deResp, ayResp, stimID, stimDir] = fetch1(theseStimInfo,...
                'stimt_df_resp','stimt_de_resp','stimt_ay_resp','stimt_targ_id','stimt_vis_dir');
            dfResp = double(dfResp); deResp = double(deResp); ayResp = double(ayResp);
            
            infDist = fetchn(v1inf.Influence & key,'inf_dist','ORDER BY targ_id, neur_id');
            infDist = reshape(infDist,size(dfResp,3),[]);

            [dfMat, dfStd] = calcAvgResps(dfResp, infDist, stimID);
            [deMat, deStd] = calcAvgResps(deResp, infDist, stimID);
            [ayMat, ayStd] = calcAvgResps(ayResp, infDist, stimID);
            
            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = nNeuron;
                    keys(nNeuron,nTarg).targ_id = nTarg;
                    keys(nNeuron,nTarg).std_df = dfStd(nNeuron,nTarg);
                    keys(nNeuron,nTarg).std_de = deStd(nNeuron,nTarg);
                    keys(nNeuron,nTarg).std_ay = ayStd(nNeuron,nTarg);
                    keys(nNeuron,nTarg).sta_df = dfMat(:,nNeuron,nTarg);
                    keys(nNeuron,nTarg).sta_de = deMat(:,nNeuron,nTarg);
                    keys(nNeuron,nTarg).sta_ay = ayMat(:,nNeuron,nTarg);
                end
            end
            
            insert(self,keys),
            
        end
    end
end

function [respMat, stdMat] = calcAvgResps(trialResps, infDist, stimID)

distThresh = 30;
preRange = 10:15; %indices for 'baseline' measurement, should be -6:-1 frames from stim onset
stdRange = 23:45; %indices for normalization, should be 7-29 frames from stim onset (exclude artifact frames)
[nRespCells,nTargets] = size(infDist);

%% Calculate Normalization Factor
for n=1:nRespCells
    validStim = find(infDist(n,:) >=distThresh);
    validInd = ismember(stimID, validStim);
    validResps = trialResps(stdRange,validInd,n);
    stdVal(n) = std(mean(validResps));
end
stdMat = repmat(stdVal(:),1,nTargets);

% Subtract baseline period for each neuron
% trialResps = trialResps-mean(trialResps(preRange,:,:));
tmp = reshape(trialResps,[],size(trialResps,3));
trialResps = trialResps - reshape(prctile(tmp,10),1,1,size(trialResps,3));

%% Calculate stim-triggered average traces
respMat = nan(size(trialResps,1),nRespCells,nTargets);
for nStimCell = 1:nTargets
    theseTri = stimID == nStimCell;
    staTrace = squeeze(mean(trialResps(:,theseTri,:),2));
    respMat(:,:,nStimCell) = staTrace;
end

end