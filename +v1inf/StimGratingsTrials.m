%{
# stimGratings Trial Projections
-> v1inf.StimGratingsData
stim_trial_id: smallint unsigned
-----
-> v1inf.Target                         # Target id for this stim trial
stim_grating_ori: tinyint unsigned     # orientation of grating drift
stim_grating_proj: blob                 # 4-element vector of true projection values
stim_grating_gain: blob                 # ... of gain-projection values
stim_grating_uni: double                # uniform projection value
%}

classdef StimGratingsTrials < dj.Computed
    	methods(Access=protected)

		function makeTuples(self, key)
            stimExpt = fetch(v1inf.StimGratingsData & key,'*');
            
            [allProj,allGain,uniProj] = computeStimProjections(stimExpt,key);
            
            validTrials = find(isfinite(stimExpt.stim_targ_id));
            nTrials = length(validTrials);
            keys = repmat(key,nTrials,1);
            for nTrial = 1:nTrials
                iTrial = validTrials(nTrial);
                keys(nTrial).stim_trial_id = nTrial;
                keys(nTrial).targ_id = stimExpt.stim_targ_id(iTrial);
                keys(nTrial).stim_grating_ori = round(mod(stimExpt.stim_vis_dir(iTrial),180));
                keys(nTrial).stim_grating_proj = allProj(iTrial,:);
                keys(nTrial).stim_grating_gain = allGain(iTrial,:);
                keys(nTrial).stim_grating_uni = uniProj(iTrial);
            end
            
            self.insert(keys),
            
		end
	end
end



function [allProj,allGain,uniProj] = computeStimProjections(stimExpt,key)

% visResid = nan(size(subResp));
%     visResid(theseTrials,:) = subResp(theseTrials,:) - visAvg(:,i)';
%     for n=1:size(distMat,1)
%         invalidStim = find(distMat(n,:) < distThresh);
%         validInd = theseTrials & ~ismember(stimID,invalidStim);
%         invalidInd = theseTrials & ismember(stimID,invalidStim);
%         visAvg(n,i) = mean(subResp(validInd,n));
%         visResid(validInd,n) = subResp(validInd,n)-visAvg(n,i);
% %         visResid(invalidInd,n) = 0;
% 
%         tmpValid = find(validInd);
%         resampledInd = tmpValid(randperm(length(tmpValid),sum(invalidInd)));
%         visResid(invalidInd,n) = subResp(resampledInd,n) - visAvg(n,i);
%     end
% end

% distThresh = 15;

%% Get Average Vectors and Residuals
deResp = stimExpt.stim_de_resp;
preResp = stimExpt.stim_de_pre;
visOri = mod(stimExpt.stim_vis_dir,180);
stimID = stimExpt.stim_targ_id;

% [nX,nY] = fetchn(v1inf.Neuron & key,'neur_xc','neur_yc','ORDER BY neur_id');
% [tX,tY,tLabel] = fetchn(v1inf.Target & key,'targ_xc','targ_yc','targ_label','ORDER BY targ_id');
tLabel = fetchn(v1inf.Target & key,'targ_label','ORDER BY targ_id');
nTargTraces = sum(tLabel>.99);
% 
% distMat = sqrt((nX-tX').^2 + (nY-tY').^2);
% excludeDist = sum(distMat<distThresh,2)>0;
excludeDist = false(size(deResp,2),1);
excludeDist(1:nTargTraces) = 1;

allOri = unique(visOri);
subResp = (deResp-preResp)./std(preResp);
subResp = subResp(:,~excludeDist);
visAvg = nan(size(subResp,2),length(allOri)); 

for i=1:length(allOri)
    theseTrials = visOri==allOri(i) & ~isnan(stimID);
    visAvg(:,i) = mean(subResp(theseTrials,:),1);
end

%% Get Projections

for i=1:size(visAvg,2)
    oriNorm(i) = norm(visAvg(:,i));
end
visNorm = visAvg./oriNorm;
uniNorm = ones(size(visAvg,1),1)/sqrt(size(visAvg,1));

onResid = nan(length(visOri),1);
uniProj = nan(length(visOri),1);
allProj = nan(length(visOri),length(allOri));
allGain = nan(length(visOri),length(allOri));
for nTrial=1:length(visOri)
    thisDirInd = find(visOri(nTrial)==allOri);
    thisProj = visNorm(:,thisDirInd);
    scaledResp = subResp(nTrial,:) ./ oriNorm(thisDirInd);
    
    uniProj(nTrial) = scaledResp * uniNorm;
    onResid(nTrial) = scaledResp * thisProj;
    allProj(nTrial,:) = scaledResp * circshift(visNorm,-thisDirInd+1,2);
    allGain(nTrial,:) = (onResid(nTrial)*thisProj') * circshift(visNorm,-thisDirInd+1,2);
 end
% uniProj = visResid * uniNorm;

end