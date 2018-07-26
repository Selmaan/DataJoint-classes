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
targ_avg_srde=NULL: float   # Average normalized residual of deconv
targ_std_srde=NULL: float   # Std of normalized residual of deconv
targ_baseline=NULL: float   # Average activity during influence blocks
%}

classdef TargetAvgCorr < dj.Imported

    methods(Access=protected)
        function makeTuples(self,key)
            validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
            validInf = (v1inf.Influence & key) & 'inf_dist>25';
            validInfRaw = (v1inf.InfVals & key) & 'inf_dist>25';
            if fetch1(v1inf.Target & key,'targ_label')>0.01
                [sD, sV] = fetchn(v1inf.AvgDeResid & (validInf & validFilt),'avg_rde','avg_sv');
                sD = sD./sV;
            else
                [sD, sV] = fetchn(v1inf.AvgLOOResid & (validInf & validFilt),'avgloo_rde','avgloo_sv');
                sD = sD./sV;
            end
            [iC, sN] = fetchn(validInf & validFilt,'inf_naivecorr','inf_shuf_n');
            sN_raw = fetchn(validInfRaw & validFilt,'inf_n_raw');
            key.targ_avg_corr = nanmean(iC);
            key.targ_std_corr = nanstd(iC);
            key.targ_avg_inf = mean(sN);
            key.targ_std_inf = std(sN);
            key.targ_avg_inf_raw = mean(sN_raw);
            key.targ_std_inf_raw = std(sN_raw);
            key.targ_avg_srde = mean(sD);
            key.targ_std_srde = std(sD);
            
            key.targ_baseline = fetchn(v1inf.AvgDeResid & (v1inf.Target & key), 'avg_mv');
            self.insert(key),
        end
    end
end