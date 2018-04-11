%{
# Data for each photostim target in an experiment
-> v1inf.Neuron
---
neur_avg_corr=NULL: float # Average Naive/Trace Correlation with other neurons
neur_std_corr=NULL: float # Std of Naive/Trace Correlation with other neurons
neur_avg_inf: float       # Average influence of this neuron
neur_std_inf: float       # Std influence of this neuron
%}

classdef NeuronAvgCorr < dj.Imported

    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)
        function makeTuples(self,key)
            validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
            validStim = proj(v1inf.Target & (v1inf.Influence & 'inf_shuf_p>5'));
            validInf = (v1inf.Influence & key) & 'inf_dist>25';
            [iC, sN, nID] = fetchn((validInf & validStim) & validFilt,...
                'inf_naivecorr','inf_shuf_n','neur_id','ORDER BY neur_id');
            allNeur = unique(nID);
            keys = repmat(key,length(allNeur),1);
            for i = 1:length(allNeur)
                keys(i).neur_id = allNeur(i);
                neurInd = nID==allNeur(i);
                keys(i).neur_avg_corr = nanmean(iC(neurInd));
                keys(i).neur_std_corr = nanstd(iC(neurInd));
                keys(i).neur_avg_inf = nanmean(sN(neurInd));
                keys(i).neur_std_inf = nanstd(sN(neurInd));
            end
            self.insert(keys(:)),
        end
    end
end