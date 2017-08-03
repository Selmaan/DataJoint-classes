%{
# Influence Measurement Directed Pair (Regression-based)
-> v1inf.Neuron
-> v1inf.Target
-----
inf_dist: double            # Pairwise distance (in um)
inf_val=NULL: double        # Influence Value
inf_naivecorr=NULL: double  # Neuron-Target Trace Correlation
%}

classdef Influence < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            theseNeurons = v1inf.Neuron & key;
            theseTargets = v1inf.Target & key;
            theseStimInfo = v1inf.StimGratingsData & key;
            theseSyncInfo = v1inf.ExpSync & key;
            
            [targ_id, targ_xc, targ_yc, neur_id] = fetchn(...
                theseTargets,'targ_id','targ_xc','targ_yc','neur_id','ORDER BY targ_id');
            [neur_id, neur_xc, neur_yc, neur_deconv] = fetchn(...
                theseNeurons,'neur_id','neur_xc','neur_yc','neur_deconv','ORDER BY neur_id');
            stimData = fetch(theseStimInfo,'*');
            syncData = fetch(theseSyncInfo,'stim_blocks','frame_times');
        end
    end
end