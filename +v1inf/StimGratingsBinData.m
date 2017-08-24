%{
# stimGratings Trials
-> v1inf.ExpSync
-> v1inf.StimGratingsBinRanges
-----
stim_targ_id: longblob
stim_vis_dir: longblob              # stim grating drift direction
stim_mv_spd: longblob               # Mouse running speed during grating presentation
stim_de_resp: longblob              # Deconvolved gratings responses for all neurons
%}

classdef StimGratingsBinData < dj.Computed
    	methods(Access=protected)

		function makeTuples(self, key)
%             q1 = sprintf('mouse_id = %d', key.mouse_id);
%             q2 = sprintf('exp_date = "%s"',key.exp_date);
            thisSync = v1inf.ExpSync & key;
            theseNeurons = v1inf.Neuron & key;
            neur_deconv = fetchn(theseNeurons,'neur_deconv','ORDER BY neur_id');
            neur_deconv = cell2mat(neur_deconv')';
            syncInfo = fetch(thisSync,'*');
            binRange = fetch1(v1inf.StimGratingsBinRanges & key,'bin_range');
            [stimID,visDir,mvSpd,deResp] = ...
                stimGratingsResp(syncInfo, neur_deconv, binRange);
            
            key.stim_targ_id = stimID;
            key.stim_vis_dir = visDir;
            key.stim_mv_spd = mvSpd;
            key.stim_de_resp = deResp;

            self.insert(key)
		end
	end
end

function [stimID,visDir,mvSpd,deResp] = ...
    stimGratingsResp(syncInfo, neur_deconv, respFrameRange)

stimID=[];visDir=[];mvSpd=[];deResp=[];

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

    visDir = cat(1,visDir,thisVisDir);
    stimID = cat(1,stimID,thisStimID);
    mvSpd = cat(1,mvSpd,thisMvSpd);
    deResp = cat(1,deResp,thisDeResp);
end
end