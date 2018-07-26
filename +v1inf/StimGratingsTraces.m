%{
# stimGratings Trials using dF/F Traces
-> v1inf.ExpSync
-----
stimt_targ_id: longblob
stimt_vis_dir: longblob              # stim grating drift direction
stimt_mv_spd: longblob               # Mouse running speed during grating presentation
stimt_de_resp: longblob              # deconvolved gratings responses for all neurons
stimt_df_resp: longblob              # dF/F gratings responses for all neurons
stimt_ay_resp: longblob              # 'raw' gratings responses for all neurons
%}

classdef StimGratingsTraces < dj.Computed
    	methods(Access=protected)

		function makeTuples(self, key)
            if strcmp(fetch1(v1inf.ExpType & key,'exp_type'), 'Multi-Contrast')
                warning('Skipping Multi-Contrast Experiment'),
                return
            end
            thisSync = v1inf.ExpSync & key;
            theseNeurons = v1inf.NeurDf & key;
            [neur_deconv, neur_df, neur_ay] = fetchn(theseNeurons,...
                'n2_deconv','n2_df','n2_raw','ORDER BY neur_id');
            neur_deconv = cell2mat(neur_deconv')';
            neur_df = cell2mat(neur_df')';
            neur_ay = cell2mat(neur_ay')';
            neur_ay = neur_ay ./ std(neur_ay,[],2) .* std(neur_df,[],2); % Normalize to scale of dF/F
            syncInfo = fetch(thisSync,'*');
            [stimID,visDir,mvSpd,deResp,dfResp,ayResp] = stimGratingsResp(...
                syncInfo, neur_deconv, neur_df, neur_ay);
            
            key.stimt_targ_id = stimID;
            key.stimt_vis_dir = visDir;
            key.stimt_mv_spd = single(mvSpd);
            key.stimt_de_resp = single(deResp);
            key.stimt_df_resp = single(dfResp);
            key.stimt_ay_resp = single(ayResp);

            self.insert(key)
		end
	end
end

function [stimID,visDir,mvSpd,deResp,dfResp,ayResp] = stimGratingsResp(...
                syncInfo, neur_deconv, neur_df, neur_ay)

stimID=[];visDir=[];mvSpd=[];deResp=[];dfResp=[];ayResp=[];

respFrameRange = -15:45;

ballVel = cat(1,syncInfo.ball_vel{:});
ballVel = bsxfun(@minus,ballVel,mode(ballVel));
for nBlock = find(syncInfo.stim_blocks)
    blockOffsetFrame = length(cat(1,syncInfo.frame_times{1:nBlock-1}));
    lastTrial = find(~isnan(syncInfo.psych2frame{nBlock}),1,'last')-1;
    visFrames = syncInfo.psych2frame{nBlock}(1:lastTrial) + blockOffsetFrame;
    thisVisDir = [syncInfo.stim_info{nBlock}((1:3:lastTrial*3)+4),...
        syncInfo.stim_info{nBlock}((2:3:lastTrial*3)+4)];
    thisStimID = syncInfo.stim_order{nBlock}';
    thisStimID(end+1:lastTrial) = nan;
    
    respFrames = bsxfun(@plus,visFrames,respFrameRange)';
    respFrames(respFrames>size(ballVel,1)) = size(ballVel,1);
    thisMvSpd = reshape(ballVel(respFrames,:),[size(respFrames), size(ballVel,2)]);
    thisMvSpd = sqrt(sum(thisMvSpd.^2,3));
    % Deconvolved
    thisDeResp = reshape(neur_deconv(:,respFrames),...
           [size(neur_deconv,1),size(respFrames)]);
    thisDeResp = permute(thisDeResp, [2 3 1]);
    % dF
    thisDfResp = reshape(neur_df(:,respFrames),...
           [size(neur_df,1),size(respFrames)]);
    thisDfResp = permute(thisDfResp, [2 3 1]);
    % AY
    thisAyResp = reshape(neur_ay(:,respFrames),...
           [size(neur_ay,1),size(respFrames)]);
    thisAyResp = permute(thisAyResp, [2 3 1]);
    
    visDir = cat(1,visDir,thisVisDir);
    stimID = cat(1,stimID,thisStimID);
    mvSpd = cat(2,mvSpd,thisMvSpd);
    deResp = cat(2,deResp,thisDeResp);
    dfResp = cat(2,dfResp,thisDfResp);
    ayResp = cat(2,ayResp,thisAyResp);
end
end