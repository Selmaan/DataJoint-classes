%{
# RandomGratings protocol trial info and responses
-> v1inf.Experiment
-----
rg_resps: longblob       # randomGratings resp structure
%}

classdef RandomGratingsExp < dj.Computed

	methods(Access=protected)

		function makeTuples(self, key)
            q1 = sprintf('mouse_id = %d', key.mouse_id);
            q2 = sprintf('exp_date = "%s"',key.exp_date);
            thisExp = v1inf.ExpSync & q1 & q2;
            
            key.rg_resps = makeRgRespStruct(thisExp);
            self.insert(key)
		end
	end

end

function allResp = makeRgRespStruct(thisExp)

allResp = struct;
allResp.nCycles = [];
dF_deconv = fetchn(v1inf.Neuron & thisExp,'neur_deconv');
dF_deconv = cat(2,dF_deconv{:})';
gratingBlocks = find(~fetch1(thisExp,'stim_blocks'));
stimInfo = fetch1(thisExp,'stim_info');
ballVel = fetch1(thisExp,'ball_vel');
psychTimes = fetch1(thisExp,'psych_times');
frameTimes = fetch1(thisExp,'frame_times');
vel_offset = mode(reshape(cat(1,ballVel{:}),[],1));
ballVel = cellfun(@(x) x-vel_offset, ballVel, 'UniformOutput', false);
for nBlock = gratingBlocks
    Dir = stimInfo{nBlock}(6:6:end);
    SF = stimInfo{nBlock}(7:6:end);
    TF = stimInfo{nBlock}(8:6:end);
    CT = stimInfo{nBlock}(9:6:end);
    
    fSpd = sqrt(sum(ballVel{nBlock}.^2,2));
    
    f2p = interp1(psychTimes{nBlock},...
        1:length(psychTimes{nBlock}),frameTimes{nBlock},'nearest');
    f2p_valid = find(~isnan(f2p));
    fCT = nan(length(frameTimes{nBlock}),1);
    fDir(f2p_valid) = Dir(f2p(f2p_valid));
    fSF(f2p_valid) = log2(SF(f2p(f2p_valid)));
    fTF(f2p_valid) = log2(TF(f2p(f2p_valid)));
    fCT(f2p_valid) = CT(f2p(f2p_valid));
    smCT = conv(fCT,gausswin(15)/sum(gausswin(15)),'same');
    ctTrig = 1+find(diff(smCT(2:end))>0 & diff(smCT(1:end-1)) < 0);
    blockOffsetFrame = length(cat(1,frameTimes{1:nBlock-1}));
    traces = dF_deconv(:,blockOffsetFrame+1:blockOffsetFrame+length(frameTimes{nBlock}));
    cycleDur = mode(diff(ctTrig));
    allResp.yPh{nBlock} = zeros(size(traces,1),cycleDur);
    fprintf('Block %d had %d cycles at %d frames-per-cycle \n',nBlock,length(ctTrig)-1,cycleDur);
    allResp.nCycles(end+1) = length(ctTrig)-1;
    for nCycle = 1:length(ctTrig)-1
        cycleInd = ctTrig(nCycle) + (1:cycleDur);
        allResp.Y{nBlock}(:,nCycle) = sum(traces(:,cycleInd),2);
        allResp.yPh{nBlock} = allResp.yPh{nBlock} + traces(:,cycleInd);
        allResp.Dir{nBlock}(nCycle) = mode(fDir(cycleInd));
        allResp.SF{nBlock}(nCycle) = mode(fSF(cycleInd));
        allResp.TF{nBlock}(nCycle) = mode(fTF(cycleInd));
        allResp.spd{nBlock}(nCycle) = mean(fSpd(cycleInd));
    end
end

allResp.Y = cat(2,allResp.Y{:})';
allResp.yPh = mean(cat(3,allResp.yPh{:}),3)';
allResp.Dir = cat(2,allResp.Dir{:})';
allResp.SF = cat(2,allResp.SF{:})';
allResp.TF = cat(2,allResp.TF{:})';
allResp.spd = cat(2,allResp.spd{:})';
       
end