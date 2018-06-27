%{
# Table with correctly scaled dF/F Traces for single-contrast experiments
-> v1inf.Neuron
---
n2_raw: longblob           # nFrames x 1 vector of non-subtracted, non-normalized trace
n2_df: longblob            # nFrames x 1 vector of cell's dF/F trace
n2_deconv: longblob        # nFrames x 1 vector of deconvolved activity
%}

classdef NeurDf < dj.Imported
    properties
        popRel = v1inf.Experiment;
    end
    
    methods(Access=protected)
        function makeTuples(self,key)
            if strcmp(fetch1(v1inf.ExpType & key,'exp_type'), 'Multi-Contrast')
                warning('Skipping Multi-Contrast Experiment'),
                return
            end
            thisDir = fetch1(v1inf.Experiment & key,'server_dir');
            expt_proc = fetch1(v1inf.Experiment & key,'expt_proc');
            
            if expt_proc == 'F'
                fprintf('Unprocessed Experiment m%d date %s \n',key.mouse_id,key.exp_date),
                return
            end
            fprintf('Loading Experiment m%d date %s \n',key.mouse_id,key.exp_date),
            if ~exist(fullfile(thisDir,'stimExpt.mat'),'file')
                thisDir = strrep(thisDir, 'Lab\','Lab\Tier1\');
                updateDir = 1;
            else
                updateDir = 0;
            end
            load(fullfile(thisDir,'stimExpt.mat')),
            load(fullfile(thisDir,'resFOV1.mat')),
            if updateDir
                resFOV1.roiInfo.slice.NMF.filename = strrep(resFOV1.roiInfo.slice.NMF.filename, 'Lab\S','Lab\Tier1\S');
                resFOV1.roiInfo.slice.NMF.traceFn = strrep(resFOV1.roiInfo.slice.NMF.traceFn, 'Lab\S','Lab\Tier1\S');
            end


            newData = deconvData(resFOV1,stimExpt);
            
            nNeurons = size(stimExpt.cellFilts,2);
            for nNeuron = 1:nNeurons
                key.neur_id = nNeuron;
                key.n2_df = newData.dF(nNeuron,:)';
                key.n2_deconv = newData.deconv(nNeuron,:)';
                key.n2_raw = newData.raw(nNeuron,:)';
                self.insert(key),
            end

            fprintf('Populated Neuron data for %d %s \n',key.mouse_id,key.exp_date);

        end
    end
end



function newData = deconvData(acqObj,stimExpt)

stimFrames = cell(0);
for nBlock = find(stimExpt.stimBlocks)
    blockOffsetFrame = length(cat(1,stimExpt.frameTimes{1:nBlock-1}));
    stimFrames{nBlock} = blockOffsetFrame + stimExpt.psych2frame{nBlock}(1:length(stimExpt.stimOrder{nBlock}));
end

stimFrames = cat(1,stimFrames{stimExpt.stimBlocks});
stimFrames = repmat(stimFrames,1,4) + repmat(0:2:6,size(stimFrames,1),1);
stimFrames = stimFrames(:);
interpFrames = cell(0);
interpFrames{1} = stimFrames;
interpFrames{2} = stimFrames+1;
interpFrames{3} = stimFrames-1;

validSources{1} = stimExpt.cIds;
[dF,deconv,denoised,Gs,Lams,A,b,f, C_raw] = extractTraces_NMF(acqObj,validSources,interpFrames);

newData.dF = cell2mat(dF);
newData.deconv = cell2mat(deconv);
newData.raw = cell2mat(C_raw);
end