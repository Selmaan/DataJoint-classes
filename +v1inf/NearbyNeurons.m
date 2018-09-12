%{
# Nearest Cell for all control & unmatched stim targets
-> v1inf.Target
-----
neur_ids: blob         # Neuron IDs nearby to this target
targ_type: enum('Matched','Unmatched','Control')    # type of photostim target
mean_dist: double       # Mean distance of selected nearby neurons
%}

classdef NearbyNeurons < dj.Computed
        
    methods(Access=protected)

        function makeTuples(self, key)
           thisTarg = fetch(v1inf.Target & key, 'targ_xc','targ_yc','targ_label','neur_id');
           [neurX,neurY,neurID] = fetchn(v1inf.Neuron & key,'neur_xc','neur_yc','neur_id');
           targDist = sqrt((neurX-thisTarg.targ_xc).^2 + (neurY-thisTarg.targ_yc).^2);
           
%            key.neur_ids = neurID(targDist<25);
           nearbyInd = targDist<=prctile(targDist,2.5);
           key.neur_ids = neurID(nearbyInd);
           key.mean_dist = mean(targDist(nearbyInd));
           if thisTarg.targ_label > 0.99
               key.targ_type = 'Matched';
           elseif thisTarg.targ_label < 0.01
               key.targ_type = 'Control';
           else
               key.targ_type = 'Unmatched';
           end
           
           self.insert(key),
            
        end
	end
end