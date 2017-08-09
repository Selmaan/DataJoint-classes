%{
# Influence Overlap Degree (processed stim image-derived)
-> v1inf.Neuron
-> v1inf.Target
-----
filt_overlap : double   # Spatial overlap of neuron filter and target stim-image
%}

classdef FiltOverlap < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            theseNeurons = v1inf.Neuron & key;
            theseTargets = v1inf.Target & key;
            [targ_procim] = fetchn(...
                theseTargets,'targ_procim','ORDER BY targ_id');
            targ_procim = cat(3,targ_procim{:});
            targ_procim(isnan(targ_procim))=0;
            [neur_filt] = fetchn(...
                theseNeurons,'neur_filt','ORDER BY neur_id');
            neur_filt = cat(2,neur_filt{:});
            neur_filt = neur_filt./sum(neur_filt);
            %%
            targ_threshim = targ_procim > 0;
            se = strel('rectangle',[3 3]);
            for i=1:size(targ_threshim,3)
                targ_threshim(:,:,i) = bwareafilt(imclose(targ_threshim(:,:,i),se),1);
            end
            targ_clipim = reshape(targ_procim.*targ_threshim,512^2,size(targ_procim,3));
            
            overlapMat = neur_filt' * targ_clipim;
            
            for nNeuron = 1:size(overlapMat,1)
                for nTarget = 1:size(overlapMat,2)
                    key.neur_id = nNeuron;
                    key.targ_id = nTarget;
                    key.filt_overlap = overlapMat(nNeuron,nTarget);
                    self.insert(key),
                end
            end
            
        end
    end
end