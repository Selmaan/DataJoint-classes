%{
# Experiment
-> v1inf.Mouse
exp_date: date                 # Date of the experiment
---
server_dir: varchar(120)         # location of data on HMS-neurobio server
expt_proc: enum('T','F')       # Whether a stimExpt file has been produced for this experiment
%}

classdef Experiment < dj.Manual
end