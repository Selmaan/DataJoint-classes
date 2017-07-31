%{
# Data for each neural cell body identified in an experiment
-> v1inf.Experiment
neur_id: smallint unsigned   # cell number
---
neur_filt: longblob              # Neural spatial filter (from NMF source extraction)
neur_xc: float               # X coordinate of neuron's centroid (in um)
neur_yc: float               # Y coordinate of neuron's centroid (in um)
neur_df: longblob            # nFrames x 1 vector of cell's dF/F trace
neur_deconv: longblob        # nFrames x 1 vector of deconvolved activity
%}

classdef Neuron < dj.Imported

    methods(Access=protected)
        function makeTuples(self,key)
            
            q1 = sprintf('mouse_id = %d', key.mouse_id);
            q2 = sprintf('exp_date = "%s"',key.exp_date);
            thisDir = fetch1(v1inf.Experiment & q1 & q2,'server_dir');
            expt_proc = fetch1(v1inf.Experiment & q1 & q2,'expt_proc');
            
            if expt_proc == 'F'
                fprintf('Unprocessed Experiment m%d date %s \n',key.mouse_id,key.exp_date),
                return
            end
            fprintf('Loading Experiment m%d date %s \n',key.mouse_id,key.exp_date),
            load(fullfile(thisDir,'stimExpt.mat')),
            
            nNeurons = size(stimExpt.cellFilts,2);
            for nNeuron = 1:nNeurons
                key.neur_id = nNeuron;
                key.neur_filt = single(full(stimExpt.cellFilts(:,nNeuron)));
                key.neur_xc = stimExpt.cellCentroids(nNeuron,1) * stimExpt.xConvFactor;
                key.neur_yc = stimExpt.cellCentroids(nNeuron,2) * stimExpt.yConvFactor;
                key.neur_df = stimExpt.dF(nNeuron,:)';
                key.neur_deconv = stimExpt.dF_deconv(nNeuron,:)';
                self.insert(key),
            end

            fprintf('Populated Neuron data for %d %s \n',key.mouse_id,key.exp_date);

        end
    end
end