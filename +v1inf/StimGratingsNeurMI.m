%{
# stimGratings Neuron Orientation Info
-> v1inf.StimGratingsData
-> v1inf.Neuron
-----
nb_ori_info: double         # mutual info between response and grating ori
nb_ori_pref: blob           # mean response to the 4 grating orientations
nb_stim_resp=NULL: double   # mean stimulation response, if neur is target
%}

classdef StimGratingsNeurMI < dj.Computed
    properties
        popRel = v1inf.StimGratingsData
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            stimExpt = fetch(v1inf.StimGratingsData & key,'*');
            theseNeur = fetch(v1inf.Neuron & key,'ORDER BY neur_id');
            nNeurons = length(theseNeur);
            if nNeurons ~= size(stimExpt.stim_de_resp,2)
                error('Neuron ID not consistent!'),
            end
            
            [neuronOriInfo,avgResp,stimResp] = computeStimMI(stimExpt,key);
        
            keys = repmat(key,nNeurons,1);
            for nNeur = 1:nNeurons
                keys(nNeur).neur_id = theseNeur(nNeur).neur_id;
                keys(nNeur).nb_ori_info = neuronOriInfo(nNeur);
                keys(nNeur).nb_ori_pref = avgResp(nNeur,:);
                keys(nNeur).nb_stim_resp = stimResp(nNeur);
            end
            
            self.insert(keys),
            
		end
	end
end



function [neuronOriInfo,avgResp,stimResp] = computeStimMI(stimExpt,key)

%% Formatting
deResp = stimExpt.stim_de_resp;
stimID = stimExpt.stim_targ_id;
% Identify zero-contrast trials, set to orientation of 255
if size(stimExpt.stim_vis_dir,2)==1
    visOri = mod(stimExpt.stim_vis_dir,180);
elseif size(stimExpt.stim_vis_dir,2)==2
    zeroContrasts = stimExpt.stim_vis_dir(:,2)==0;
    if ~sum(zeroContrasts)
        error('Data Formatting error, see above in code?'),
    end
    visOri = mod(stimExpt.stim_vis_dir(:,1),180);
    visOri(zeroContrasts) = 255;
else
    error('Unknown Data Format'),
end

[tID,tLabel] = fetchn(v1inf.Target & key,'targ_id','targ_label','ORDER BY targ_id');
neurID2targID = tID(tLabel>.99);
allOri = unique(visOri);

%% Bin Responses

binResp = nan(size(deResp));
avgResp = nan(size(deResp,2),length(allOri));
stimResp = nan(size(deResp,2),1);
for nNeur = 1:size(binResp,2)
    if nNeur<=length(neurID2targID)
        thisTarg = neurID2targID(nNeur);
        validTrials = find(stimID~=thisTarg);
        stimResp(nNeur) = mean(deResp(stimID==thisTarg,nNeur));
    else
        validTrials = 1:length(stimID);
    end
    y = deResp(validTrials,nNeur);
    
    for iOri = 1:length(allOri)
        avgResp(nNeur,iOri) = mean(y(visOri(validTrials)==allOri(iOri)));
    end
    
    binResp(validTrials(y==0),nNeur) = 0;
    yBins = [0, prctile(y(y>0),[25 50 75]), inf];
    for nBin = 1:length(yBins)-1
        ind = y>yBins(nBin) & y<=yBins(nBin+1);
        binResp(validTrials(ind),nNeur) = nBin;
    end
end
binIDs = unique(binResp(~isnan(binResp(:))));

%% Calculate Base and Condition Probabilities and Neuron-Ori Mutual-Info
baseProb = nan(size(binResp,2),length(binIDs));
condProb = nan(size(binResp,2),length(binIDs),length(allOri));
for nBin = 1:length(binIDs)
    iBin = binIDs(nBin);
    baseProb(:,nBin) = sum(binResp==iBin)./sum(~isnan(binResp));
    for nOri = 1:length(allOri)
        ind = visOri==allOri(nOri);
        condProb(:,nBin,nOri) = sum(binResp(ind,:)==iBin)./sum(~isnan(binResp(ind,:)));
    end
end
condProb(condProb==0) = 1/(1+size(binResp,1));

% sum ratio of probability of conditional versus base prob over all bins of response and orientation
neuronOriInfo = sum(sum(condProb.*log2(condProb./baseProb),3),2);


end