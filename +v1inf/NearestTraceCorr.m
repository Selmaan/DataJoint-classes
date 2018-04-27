%{
# NearestNeuron Trace Correlations
-> v1inf.Neuron
-> v1inf.TargetedNeurons
-----
nearest_dist: double            # Pairwise distance (in um)
nearest_tracecorr: double  # Neuron-Target Trace Correlation
%}

classdef NearestTraceCorr < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            theseNeurons = v1inf.Neuron & key;
            theseTargets = v1inf.TargetedNeurons & key;
            theseSyncInfo = v1inf.ExpSync & key;
            
              % Get locations from neuron table
            [targ_id, targ_xc, targ_yc, targ_neur_id] = fetchn(theseTargets*v1inf.Neuron,...
                'targ_id','neur_xc','neur_yc','neur_id','ORDER BY targ_id');
        
            [neur_xc, neur_yc, neur_deconv] = fetchn(...
                theseNeurons,'neur_xc','neur_yc','neur_deconv','ORDER BY neur_id');
            syncData = fetch(theseSyncInfo,'stim_blocks','frame_times');
            
            infDist = sqrt((neur_xc-targ_xc').^2 + (neur_yc-targ_yc').^2);
            naiveCorrs = computeNaiveCorr(neur_deconv, targ_neur_id, syncData);
            naiveCorrs(~isfinite(naiveCorrs))=0;
                        
            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = nNeuron;
                    keys(nNeuron,nTarg).targ_id = nTarg;
                    keys(nNeuron,nTarg).nearest_dist = infDist(nNeuron,nTarg);
                    keys(nNeuron,nTarg).nearest_tracecorr = naiveCorrs(nNeuron,nTarg);
                end
            end
            
            insert(self,keys),
            
        end
    end
end

function naiveCorrs = computeNaiveCorr(neur_deconv, targ_neur_id, syncData)


neur_deconv = cell2mat(neur_deconv');
nNeurons = size(neur_deconv,2);
nTargets = length(targ_neur_id);
validTarg = find(~isnan(targ_neur_id));

validBlocks = find(~syncData.stim_blocks);
naiveCorrs = nan(nNeurons,nTargets,length(validBlocks));
for nBlock = validBlocks
    blockOffsetFrame = length(cat(1,syncData.frame_times{1:nBlock-1}));
    blockInd = blockOffsetFrame + (1:length(syncData.frame_times{nBlock}));
    bin_deconv = squeeze(mean(reshape(neur_deconv(blockInd,:),10,[],nNeurons),1));
    naiveCorrs(:,validTarg,validBlocks==nBlock) = corr(bin_deconv,bin_deconv(:,targ_neur_id(validTarg)));
end
naiveCorrs = nanmean(naiveCorrs,3);
end