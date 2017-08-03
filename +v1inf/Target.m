%{
# Data for each photostim target in an experiment
-> v1inf.Experiment
targ_id: smallint unsigned      # photostim target number
---
targ_xc: float                  # target X coordinate (in um)
targ_yc: float                  # target Y coordinate (in um)
targ_rawim: longblob            # raw dF/F stim response image
targ_procim: longblob           # processed stim response image
targ_label: float               # target label (0=control, 0.1=no ID, 1=good, 1.2=messy or partially split/merged source)
-> [nullable] v1inf.Neuron      # Target neuron/trace id, if it exists
%}

classdef Target < dj.Imported

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
            
            nTargets = size(stimExpt.roiCentroid,1);
            for nTarget = 1:nTargets
                key.targ_id = nTarget;
                key.targ_xc = stimExpt.roiCentroid(nTarget,1) * stimExpt.xConvFactor;
                key.targ_yc = stimExpt.roiCentroid(nTarget,2) * stimExpt.yConvFactor;
                key.targ_rawim = single(stimExpt.rawStimIm(:,:,nTarget));
                key.targ_procim = single(stimExpt.procStimIm(:,:,nTarget));
                key.targ_label = stimExpt.targetLabel(nTarget);
                if key.targ_label >=1
                   key.neur_id = sum(stimExpt.targetLabel(1:nTarget)>=1);
                else
                   key.neur_id = [];
                end
                self.insert(key),
            end

            fprintf('Populated Target data for %d %s \n',key.mouse_id,key.exp_date);

        end
    end
end