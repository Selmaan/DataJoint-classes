%{
# ExpType
-> v1inf.Experiment
---
exp_type: enum('Single-Contrast','Multi-Contrast','Monitor-Off','Single-Contrast-Dense') # Experiment protocol and expression density
%}

classdef ExpType < dj.Manual
end


%% Experiment Type listing
% Experiments w/ mice 18-29 are all single-contrast (except for two in m26)
% Experiments w/ mouse 33 are single-contrast-dense 
% Experiments w/ mouse 37 from 12-13 to 12-27  are monitor-off (does not include those in 2018)
% Experiments w/ mouse 37 in 2018 and mice 38-41 are multi-contrast