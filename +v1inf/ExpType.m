%{
# ExpType
-> v1inf.Experiment
---
exp_type: enum('Single-Contrast','Multi-Contrast','Monitor-Off','Single-Contrast-Dense') # Experiment protocol and expression density
%}

classdef ExpType < dj.Manual
end


%% Experiment Type listing
% Experiments w/ mice 18-29 are all single-contrast
% Experiments w/ mouse 33 are single-contrast-dense 
% Experiments w/ mouse 37 from 12-13 to 12-27  are monitor-off (does not include those in Jan. 2018)