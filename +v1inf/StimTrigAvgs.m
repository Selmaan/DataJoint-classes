%{
# stimGratings Trials
-> v1inf.Target
-----
sta_image: longblob         # stim-triggered-average image
sta_mov_fn: varchar(120)    # filename of stim-triggered-average movie
%}

classdef StimTrigAvgs < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
        
    methods(Access=protected)

        function makeTuples(self, key)
            % Set parameters
            offsetFrames = setdiff(-24:60,[0:2:6, 31:2:37]);
            postStimInd = 28:51; %7-30 frames after stimulation
            preStimInd = 1:24; %1-24 frames before stimulation
            
            %TEMPORARY% - populate monitor-off experiments only
%             thisExptType = fetch1(v1inf.ExpType & key,'exp_type');
%             if ~strcmp(thisExptType,'Monitor-Off')
%                 return
%             end
            
            % Get Data
            syncInfo = fetch(v1inf.ExpSync & key, ...
                'stim_blocks','frame_times','psych2frame','stim_order');
            fn = fetch1(v1inf.Experiment & key,'server_dir');
            newFn = false;
            if ~exist(fn,'dir')
                fprintf('Updating File Locations to Tier 1 \n'),
                fn = strrep(fn,'\Selmaan','\Tier1\Selmaan');
                newFn = true;
            end
            load(fullfile(fn,'resFOV1.mat')),
            % Update corrected movie locations if needed
            if newFn
                for i = 1:length(resFOV1.correctedMovies.slice.channel.fileName)
                    oldFn = resFOV1.correctedMovies.slice.channel.fileName{i};
                    newFn = strrep(oldFn,'\Selmaan','\Tier1\Selmaan');
                    resFOV1.correctedMovies.slice.channel.fileName{i} = newFn;
                end
            end
            
            % Collect Stimulation times for each target
            nTargs = unique([syncInfo.stim_order{:}]);
            targStimFrames = [];
            for nBlock = find(syncInfo.stim_blocks)
                blockOffsetFrame = length(cat(1,syncInfo.frame_times{1:nBlock-1}));
                thisStimID = syncInfo.stim_order{nBlock}';
                visFrames = syncInfo.psych2frame{nBlock}(1:length(thisStimID)) + blockOffsetFrame;
                theseStimFrames = [];
                for nTarg = nTargs
                    theseStimFrames(:,nTarg) = visFrames(thisStimID==nTarg);
                end
                targStimFrames = cat(1,targStimFrames,theseStimFrames);
            end

            % Create average-bin cell array
            allFrames =  cell(length(offsetFrames),length(nTargs));
            for nTarg = nTargs
                for i=1:length(offsetFrames)
                    allFrames{i,nTarg} = targStimFrames(:,nTarg) + offsetFrames(i);
                end
            end
            
            [rawMov,dFmov] = eventTriggeredMovie(resFOV1,reshape(allFrames,[],1));
            rawMov = reshape(rawMov,512,512,length(offsetFrames),length(nTargs));
            dFmov = reshape(dFmov,512,512,length(offsetFrames),length(nTargs));
            subIms = squeeze(mean(dFmov(:,:,postStimInd,:),3) - ...
                mean(dFmov(:,:,preStimInd,:),3));
            
            mkdir(fn);
            for nTarg = nTargs
                thisFile = sprintf('Targ %d STA.tif',nTarg);
                key.sta_image = subIms(:,:,nTarg);
                key.sta_mov_fn = fullfile(fn,thisFile);
                key.targ_id = nTarg;
                tiffWrite(rawMov(:,:,:,nTarg),thisFile,...
                    fn,'int16');
%                 thisDfMov = dFmov(:,:,:,nTarg);
%                 save(thisFn,'thisRawMov','thisDfMov'),
                self.insert(key),
            end
        end
	end
end