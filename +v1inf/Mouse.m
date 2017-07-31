%{
# mouse
mouse_id: smallint unsigned     # unique mouse id
---
doi: date                      # date of viral injection
inj_chr: enum('C1V1', 'ChrR')   # injected channelrhodopsin - C1V1 or ChrimsonR
%}

classdef Mouse < dj.Manual
end