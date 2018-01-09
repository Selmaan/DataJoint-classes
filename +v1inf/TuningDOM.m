%{
# Pre-calculated tuning depth-of-modulation values
-> v1inf.TuningProps
-----
dir_dom=NULL : double       # Depth-of-modulation of direction/orientation tuning
sf_dom=NULL : double        # DOM for SF
tf_dom=NULL : double        # DOM for TF
spd_dom=NULL : double       # DOM for Running Speed
%}

classdef TuningDOM < dj.Computed

	methods(Access=protected)

		function makeTuples(self, key)
            
            % Get Neuron Tuning
            [dirMean,sfMean,tfMean,spdMean,dirStd,sfStd,tfStd,spdStd] = ...
                fetch1(v1inf.TuningProps & key,'dir_mean','sf_mean','tf_mean','spd_mean',...
                'dir_std','sf_std','tf_std','spd_std');
            
            % Get DOMs for each tuning type
            key.dir_dom = tuningDOM(dirMean,dirStd);
            key.sf_dom = tuningDOM(sfMean,sfStd);
            key.tf_dom = tuningDOM(tfMean,tfStd);
            key.spd_dom = tuningDOM(spdMean,spdStd);
            
            
            self.insert(key)
		end
	end

end

function DOM = tuningDOM(tuneMean,tuneStd)

% get max and min index, using std offset to avoid choosing an overly noisey extremum
[~,maxInd] = max(tuneMean-tuneStd);
[~,minInd] = min(tuneMean+tuneStd);

% Calculate DOM
maxVal = tuneMean(maxInd);
minVal = tuneMean(minInd);
DOM = (maxVal-minVal)/sqrt(tuneStd(maxInd)^2 + tuneStd(minInd)^2);

end