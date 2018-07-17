%{
# RandomGratings protocol response magnitudes
-> v1inf.Neuron
-----
rg_mag: longblob       # randomGratings responses
%}

classdef RandomGratingsRespMag < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
	methods(Access=protected)

		function makeTuples(self, key)
            thisExpt = fetch1(v1inf.ExpType & key,'exp_type');
            if ~contains(thisExpt,'Single-Contrast')
                fprintf('Skipping all but Single-Contrast Experiments \n');
                return
            end
            thisExp = v1inf.ExpSync & key;
            dF = fetchn(v1inf.NeurDf & thisExp,'n2_df','ORDER BY neur_id');
            dF = cat(2,dF{:})';
            nNeurs = size(dF,1);
            
            Y = makeRgRespStruct(thisExp,dF);
            keys = repmat(key,nNeurs,1);
            for nNeur = 1:nNeurs
                keys(nNeur).neur_id = nNeur;
                keys(nNeur).rg_mag = Y(:,nNeur);
            end
            self.insert(keys)
		end
	end

end

function Y = makeRgRespStruct(thisExp,dF)

gratingBlocks = find(~fetch1(thisExp,'stim_blocks'));
stimInfo = fetch1(thisExp,'stim_info');
psychTimes = fetch1(thisExp,'psych_times');
frameTimes = fetch1(thisExp,'frame_times');
for nBlock = gratingBlocks
    CT = stimInfo{nBlock}(9:6:end);
        
    f2p = interp1(psychTimes{nBlock},...
        1:length(psychTimes{nBlock}),frameTimes{nBlock},'nearest');
    f2p_valid = find(~isnan(f2p));
    fCT = nan(length(frameTimes{nBlock}),1);
    fCT(f2p_valid) = CT(f2p(f2p_valid));
    smCT = conv(fCT,gausswin(15)/sum(gausswin(15)),'same');
    ctTrig = 1+find(diff(smCT(2:end))>0 & diff(smCT(1:end-1)) < 0);
    blockOffsetFrame = length(cat(1,frameTimes{1:nBlock-1}));
    
    traces = dF(:,blockOffsetFrame+1:blockOffsetFrame+length(frameTimes{nBlock}));
    traces_smoothed = sgolayfilt(traces,5,61,[],2);
    cycleDur = mode(diff(ctTrig));
    fprintf('Block %d had %d cycles at %d frames-per-cycle \n',nBlock,length(ctTrig)-1,cycleDur);
    for nCycle = 1:length(ctTrig)-1
        cycleInd = ctTrig(nCycle) + (1:cycleDur);
        Y{nBlock}(:,nCycle) = prctile(traces_smoothed(:,cycleInd),99,2)-...
            prctile(traces_smoothed(:,cycleInd),1,2);
    end
end

Y = cat(2,Y{:})';
       
end