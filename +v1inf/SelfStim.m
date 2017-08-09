%{
# Self-Stimulation Magnitude (influence regression derived)
-> v1inf.Target
-----
self_stim=NULL : double   # dF/F stim response magnitude estimate ('rawmu' field from Influence) 
%}

classdef SelfStim < dj.Computed
    methods(Access=protected)

		function makeTuples(self, key)
            thisTarget = fetch(v1inf.Target & key,'neur_id');
            if ~isnan(thisTarget.neur_id)
                key.self_stim = fetch1(v1inf.Influence & thisTarget,'inf_rawmu');
            else
                key.self_stim = nan;
            end
            
            self.insert(key),
        end
    end
end