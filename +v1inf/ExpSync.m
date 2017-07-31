%{
# Synchronization Info for Imaging, Stimulation and PsychToolbox
-> v1inf.Experiment
---
stim_blocks: tinyblob   # logical vector indicating number and type of blocks
frame_times: longblob   # 1xnBlock cell array with frame times from pClamp
psych_times: longblob   # 1xnBlock cell array with psych times from pClamp
psych2frame: longblob  # 1xnBlock cell array with psych2frame mapping
stim_order: longblob    # 1xnBlock cell array with optoStim targets for each stim block (empty for gratings blocks)
stim_info: longblob     # 1xnBlock cell array with psychtoolbox grating info
ball_vel: longblob      # 1xnBlock cell array with ball velocity for each imaging frame
%}

classdef ExpSync < dj.Imported

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
            
            key.stim_blocks = stimExpt.stimBlocks;
            key.frame_times = stimExpt.frameTimes;
            key.psych_times = stimExpt.psychTimes;
            key.psych2frame = stimExpt.psych2frame;
            key.stim_order = stimExpt.stimOrder;
            key.ball_vel = stimExpt.ballVel;
            key.stim_info = stimExpt.stimInfo;

            % insert the key into self
            self.insert(key)

            fprintf('Populated sync data for %d %s \n',key.mouse_id,key.exp_date);

        end
    end
end