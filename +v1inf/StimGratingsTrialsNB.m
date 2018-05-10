%{
# stimGratings Single Trial Naive Bayes Decoder
-> v1inf.StimGratingsData
stim_trial_id: smallint unsigned
-----
-> v1inf.Target                     # Target id for this stim trial
nb_grating_ori: tinyint unsigned    # orientation of grating drift
nb_grating_prob_r: blob               # 4-element vector of decoder log-probabilities (LOO-repetition)
nb_grating_pred_r: tinyint unsigned   # Decoder predicted orientation (LOO-repetition)
nb_grating_prob_t: blob               # 4-element vector of decoder log-probabilities (LOO-target)
nb_grating_pred_t: tinyint unsigned   # Decoder predicted orientation (LOO-target)
%}

classdef StimGratingsTrialsNB < dj.Computed
    	methods(Access=protected)

		function makeTuples(self, key)
            stimExpt = fetch(v1inf.StimGratingsData & key,'*');
            
            [cvProbs_r,cvPred_r,cvProbs_t,cvPred_t] = computeStimProbs(stimExpt,key);
            
            validTrials = find(isfinite(stimExpt.stim_targ_id));
            if length(validTrials)~=length(cvPred_r)
                error('Trial indexes mixed up!'),
            end
            
            nTrials = length(validTrials);
            keys = repmat(key,nTrials,1);
            for nTrial = 1:nTrials
                iTrial = validTrials(nTrial);
                keys(nTrial).stim_trial_id = nTrial;
                keys(nTrial).targ_id = stimExpt.stim_targ_id(iTrial);
                keys(nTrial).nb_grating_ori = round(mod(stimExpt.stim_vis_dir(iTrial),180));
                keys(nTrial).nb_grating_prob_r = cvProbs_r(nTrial,:);
                keys(nTrial).nb_grating_pred_r = cvPred_r(nTrial,:);
                keys(nTrial).nb_grating_prob_t = cvProbs_t(nTrial,:);
                keys(nTrial).nb_grating_pred_t = cvPred_t(nTrial,:);
            end
            
            self.insert(keys),
            
		end
	end
end



function [cvProbs_r,cvPred_r,cvProbs_t,cvPred_t] = computeStimProbs(stimExpt,key)

%% Formatting
deResp = stimExpt.stim_de_resp;
stimID = stimExpt.stim_targ_id;
validTrials = ~isnan(stimID);
% Identify zero-contrast trials, set to orientation of 255
if size(stimExpt.stim_vis_dir,2)==1
    visOri = mod(stimExpt.stim_vis_dir,180);
elseif size(stimExpt.stim_vis_dir,2)==2
    zeroContrasts = stimData.stim_vis_dir(:,2)==0;
    if ~sum(zeroContrasts)
        error('Data Formatting error, see above in code?'),
    end
    visOri = mod(stimExpt.stim_vis_dir,180);
    visOri(zeroContrasts) = 255;
else
    error('Unknown Data Format'),
end

tLabel = fetchn(v1inf.Target & key,'targ_label','ORDER BY targ_id');
nTargTraces = sum(tLabel>.99);
excludeDist = false(size(deResp,2),1);
excludeDist(1:nTargTraces) = 1;
allOri = unique(visOri);

validResp = deResp(validTrials,~excludeDist);
ori = visOri(validTrials);
validStim = stimID(validTrials)';
validReps = ceil((1:length(validStim))/length(unique(validStim)));
%% Bin Responses

binResp = nan(size(validResp));
for nNeur = 1:size(binResp,2)
    y = validResp(:,nNeur);
    binResp(y==0,nNeur) = 0;
    yBins = [0, prctile(y(y>0),[25 50 75]), inf];
    for nBin = 1:length(yBins)-1
        ind = y>yBins(nBin) & y<=yBins(nBin+1);
        binResp(ind,nNeur) = nBin;
    end
end
binIDs = unique(binResp(:));

%% Build LOO-Repetition Decoder

cvProbs_r = nan(sum(validTrials),length(allOri)+1);
cvPred_r = nan(sum(validTrials),1);
for nRep = unique(validReps)
    % Get Training and Testing Splits
    trainInd = validReps~=nRep;
    trainResp = binResp(trainInd,:);
    trainOri = ori(trainInd);
    testInd = validReps==nRep;
    testResp = binResp(testInd,:);

    % Get Conditional Probabilities
    baseProb = nan(size(binResp,2),length(binIDs));
    condProb = nan(size(binResp,2),length(binIDs),length(allOri));
    for nBin = 1:length(binIDs)
        iBin = binIDs(nBin);
        baseProb(:,nBin) = mean(trainResp==iBin);
        for nOri = 1:length(allOri)
            ind = trainOri==allOri(nOri);
            condProb(:,nBin,nOri) = mean(trainResp(ind,:)==iBin);
        end
    end
    condProb(condProb==0) = 1/(1+size(trainResp,1));
    
    % Get test predictions and log-likelihoods
    trialProbs = nan([size(testResp), length(allOri)+1]);
    for nNeuron = 1:size(testResp,2)
        trialProbs(:,nNeuron,:) = cat(3,...
            condProb(nNeuron,testResp(:,nNeuron)+1,:),baseProb(nNeuron,testResp(:,nNeuron)+1));
    end
    logSumProb = squeeze(sum(log(trialProbs),2));
    [~,trialPred] = max(logSumProb(:,1:length(allOri)),[],2);
    oriPred = allOri(trialPred);
    
    % Add test values to master arrays
    cvProbs_r(testInd,:) = logSumProb;
    cvPred_r(testInd,:) = oriPred;
end


%% Build LOO-Target Decoder

cvProbs_t = nan(sum(validTrials),length(allOri)+1);
cvPred_t = nan(sum(validTrials),1);
for nTarg = unique(validStim)
    % Get Training and Testing Splits
    trainInd = validStim~=nTarg;
    trainResp = binResp(trainInd,:);
    trainOri = ori(trainInd);
    testInd = validStim==nTarg;
    testResp = binResp(testInd,:);

    % Get Conditional Probabilities
    baseProb = nan(size(binResp,2),length(binIDs));
    condProb = nan(size(binResp,2),length(binIDs),length(allOri));
    for nBin = 1:length(binIDs)
        iBin = binIDs(nBin);
        baseProb(:,nBin) = mean(trainResp==iBin);
        for nOri = 1:length(allOri)
            ind = trainOri==allOri(nOri);
            condProb(:,nBin,nOri) = mean(trainResp(ind,:)==iBin);
        end
    end
    condProb(condProb==0) = 1/(1+size(trainResp,1));
    
    % Get test predictions and log-likelihoods
    trialProbs = nan([size(testResp), length(allOri)+1]);
    for nNeuron = 1:size(testResp,2)
        trialProbs(:,nNeuron,:) = cat(3,...
            condProb(nNeuron,testResp(:,nNeuron)+1,:),baseProb(nNeuron,testResp(:,nNeuron)+1));
    end
    logSumProb = squeeze(sum(log(trialProbs),2));
    [~,trialPred] = max(logSumProb(:,1:length(allOri)),[],2);
    oriPred = allOri(trialPred);
    
    % Add test values to master arrays
    cvProbs_t(testInd,:) = logSumProb;
    cvPred_t(testInd,:) = oriPred;
end


%% Cycle trial probabilities so that 1st entry is presented grating orientation

for nTrial = 1:sum(validTrials)
    oriInd = find(ori(nTrial)==allOri);
    cvProbs_r(nTrial,1:4) = circshift(cvProbs_r(nTrial,1:4),-oriInd+1,2);
    cvProbs_t(nTrial,1:4) = circshift(cvProbs_t(nTrial,1:4),-oriInd+1,2);
end

end