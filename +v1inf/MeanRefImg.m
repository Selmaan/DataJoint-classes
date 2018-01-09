%{
# Mean Reference Images
-> v1inf.Experiment
---
mean_ref: longblob # Mean Reference Image
%}

classdef MeanRefImg < dj.Imported

    methods(Access=protected)
        function makeTuples(self,key)
                        
            tmpDir = fetch1(v1inf.Experiment & key,'server_dir');
            thisDir = strrep(tmpDir,'HarveyLab\Selmaan','HarveyLab\Tier1\Selmaan');
            fprintf('\n Loading Acq Obj...'),
            load(fullfile(thisDir,'resFOV1.mat')),
            fprintf('Done!'),
            
            key.mean_ref = meanRef(resFOV1);

            % insert the key into self
            self.insert(key)

        end
    end
end