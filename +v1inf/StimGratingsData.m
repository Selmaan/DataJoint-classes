%{
# stimGratings Trials
-> v1inf.ExpSync
-----
stim_targ_id: longblob
stim_vis_dir: longblob              # stim grating drift direction
stim_mv_spd: longblob               # Mouse running speed during grating presentation
stim_pre_spd: longblob              # Mouse running speed before grating
stim_de_resp: longblob              # Deconvolved gratings responses for all neurons
stim_de_pre: longblob               # Deconvolved pre-gratings response for all neurons
%}

classdef StimGratingsData < dj.Computed
    	methods(Access=protected)

		function makeTuples(self, key)
            q1 = sprintf('mouse_id = %d', key.mouse_id);
            q2 = sprintf('exp_date = "%s"',key.exp_date);
            thisSync = v1inf.ExpSync & q1 & q2;
            theseNeurons = v1inf.Neuron & q1 & q2;
            neur_deconv = fetchn(theseNeurons,'neur_deconv','ORDER BY neur_id');
            neur_deconv = cell2mat(neur_deconv')';
            syncInfo = fetch(thisSync,'*');
            [stimID,visDir,mvSpd,mvPre,deResp,dePre] = stimGratingsResp(syncInfo, neur_deconv);
            
            key.stim_targ_id = stimID;
            key.stim_vis_dir = visDir;
            key.stim_mv_spd = mvSpd;
            key.stim_pre_spd = mvPre;
            key.stim_de_resp = deResp;
            key.stim_de_pre = dePre;

            self.insert(key)
		end
	end
end

function [stimID,visDir,mvSpd,mvPre,deResp,dePre] = stimGratingsResp(syncInfo, neur_deconv)

stimID=[];visDir=[];mvSpd=[];mvPre=[];deResp=[];dePre=[];

respFrameRange = 0:10;
preFrameRange = -14:-4;

ballVel = cat(1,syncInfo.ball_vel{:});
ballVel = bsxfun(@minus,ballVel,mode(ballVel));
for nBlock = find(syncInfo.stim_blocks)
    blockOffsetFrame = length(cat(1,syncInfo.frame_times{1:nBlock-1}));
    lastTrial = find(~isnan(syncInfo.psych2frame{nBlock}),1,'last')-1;
    visFrames = syncInfo.psych2frame{nBlock}(1:lastTrial) + blockOffsetFrame;
    thisVisDir = syncInfo.stim_info{nBlock}((1:3:lastTrial*3)+4);
    thisStimID = syncInfo.stim_order{nBlock}';
    thisStimID(end+1:length(thisVisDir)) = nan;
    
    respFrames = bsxfun(@plus,visFrames,respFrameRange)';
    thisMvSpd = reshape(ballVel(respFrames,:),[size(respFrames), size(ballVel,2)]);
    thisMvSpd = sqrt(sum(squeeze(mean(thisMvSpd,1)).^2,2));
    thisDeResp = reshape(neur_deconv(:,respFrames),...
           [size(neur_deconv,1),size(respFrames)]);
    thisDeResp = squeeze(mean(thisDeResp,2))';
    
    preFrames = bsxfun(@plus,visFrames,preFrameRange)';
    thisPreResp = reshape(neur_deconv(:,preFrames),...
           [size(neur_deconv,1),size(preFrames)]);
    thisPreResp = squeeze(mean(thisPreResp,2))';
    thisPreMvSpd = reshape(ballVel(preFrames,:),[size(preFrames), size(ballVel,2)]);
    thisPreMvSpd = sqrt(sum(squeeze(mean(thisPreMvSpd,1)).^2,2));
    
    visDir = cat(1,visDir,thisVisDir);
    stimID = cat(1,stimID,thisStimID);
    mvSpd = cat(1,mvSpd,thisMvSpd);
    mvPre = cat(1,mvPre,thisPreMvSpd);
    deResp = cat(1,deResp,thisDeResp);
    dePre = cat(1,dePre,thisPreResp);
end
end