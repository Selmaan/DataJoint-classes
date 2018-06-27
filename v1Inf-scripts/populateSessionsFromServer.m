% This script helps add new experimental sessions or mice to the database

%% Add new Sessions

% mouseIDs = fetchn(v1inf.Mouse,'mouse_id');
mouseIDs = 38:41; %Insert new mouse ID here!

for nID = 1:length(mouseIDs)
    id = mouseIDs(nID);
    baseDir = sprintf('%s%d%s','Z:\HarveyLab\Tier1\Selmaan\V1-RF\m',id);
    seshFolders = dir(sprintf('%s%s',baseDir,'\1*'));
    for nFolder = 1:length(seshFolders)
        thisDir = fullfile(baseDir,seshFolders(nFolder).name);
        if length(dir(fullfile(thisDir,'stimExpt.mat')))==1
            thisStimProc = 'T';
        elseif length(dir(fullfile(thisDir,'stimExpt.mat')))==0
            thisStimProc = 'F';
        else
            warning('Invalid Folder or File structure?'),
            keyboard,
        end
        
        tmpDate = seshFolders(nFolder).name;
        thisDate = sprintf('20%s-%s-%s',tmpDate(1:2),tmpDate(3:4),tmpDate(5:6));
        tmpPrompt = sprintf('Is experiment from m%d on %s valid?',id,tmpDate);
        includeSession = input(tmpPrompt);
        if includeSession
            insert(v1inf.Experiment,{id,thisDate,thisDir,thisStimProc}),
        end
    end
end
  
%% Specify Experiment 'Type'
newExpType = 'Multi-Contrast';
newSessions = v1inf.Experiment - v1inf.ExpType;
newKeyes = fetch(newSessions);
for i=1:length(newKeyes)
    thisKey = newKeyes(i);
    thisKey.exp_type = newExpType;
    insert(v1inf.ExpType,thisKey),
end

%%
% After adding new session to database from sever, run following code
populate(v1inf.ExpSync),
populate(v1inf.MeanRefImg),
populate(v1inf.Neuron),
populate(v1inf.Target),
populate(v1inf.TargetedNeurons),
populate(v1inf.RandomGratingsExp),
% Manually export randomGratings data for fitting on orchestra
populate(v1inf.StimGratingsData),
populate(v1inf.FiltOverlap),
populate(v1inf.Influence),
populate(v1inf.InfVals),
populate(v1inf.NearestTraceCorr),
populate(v1inf.NeuronAvgCorr),
populate(v1inf.TargetAvgCorr),
populate(v1inf.StimGratingsTrials),
populate(v1inf.StimGratingsTrialsNB),
populate(v1inf.StimGratingsNeurMI),

%%
% After fitting randomGratings on Orchestra, execute the following code
% Need to manually set up the file locations to load in GP models...
populate(v1inf.RandomGratingsGP),
populate(v1inf.TuningProps),
populate(v1inf.PairCorrs),
populate(v1inf.NearestPairCorrs),
populate(v1inf.TuningDOM),

%%
% Process STAs on their own timescale, can take quite a while...
populate(v1inf.StimTrigAvgs),

%% Time-bin Analysis also can take quite a while, tho faster than the STAs
populate(v1inf.StimGratingsBinData),
populate(v1inf.InfluenceBin),

%% Analyses & Tables for paper review process, excluding multi-contrast!
populate(v1inf.NeurDf),
populate(v1inf.StimGratingsTraces),
populate(v1inf.InfTrace),
populate(v1inf.InfTraceResid)