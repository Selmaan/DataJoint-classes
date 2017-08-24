%{
# stimGratings Response Epoch
bin_id: smallint unsigned     # bin id
---
bin_range: blob                      # frame range post-stim for binning responses
%}

classdef StimGratingsBinRanges < dj.Lookup
end