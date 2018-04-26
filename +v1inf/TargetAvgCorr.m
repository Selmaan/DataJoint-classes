%{
# Data for each photostim target in an experiment
-> v1inf.Target
---
targ_avg_corr=NULL: float # Average Naive/Trace Correlation with other neurons
targ_std_corr=NULL: float # Std of Naive/Trace Correlation with other neurons
targ_avg_inf: float       # Average influence of this target (for all valid neurons)
targ_std_inf: float       # Std of influence of this target (for all valid neurons)
targ_avg_inf_raw: float   # Average influence using 'uncentered' inf vals
targ_std_inf_raw: float   # Std of influence using 'uncentered' inf vals
%}

classdef TargetAvgCorr < dj.Imported

    methods(Access=protected)
        function makeTuples(self,key)
            validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
            validInf = (v1inf.Influence & key) & 'inf_dist>25';
            validInfRaw = (v1inf.InfVals & key) & 'inf_dist>25';
            [iC, sN] = fetchn(validInf & validFilt,'inf_naivecorr','inf_shuf_n');
            sN_raw = fetchn(validInfRaw & validFilt,'inf_n_raw');
            key.targ_avg_corr = nanmean(iC);
            key.targ_std_corr = nanstd(iC);
            key.targ_avg_inf = nanmean(sN);
            key.targ_std_inf = nanstd(sN);
            key.targ_avg_inf_raw = nanmean(sN_raw);
            key.targ_std_inf_raw = nanstd(sN_raw);
            self.insert(key),
        end
    end
end