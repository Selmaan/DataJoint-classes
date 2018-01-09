% This script helps add new experimental sessions or mice to the database

%% Add new Sessions

% mouseIDs = fetchn(v1inf.Mouse,'mouse_id');
mouseIDs = 37; %Insert new mouse ID here!

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
        if str2num(tmpDate(1:2))==17 & str2num(tmpDate(3:4))>2
            tmpPrompt = sprintf('Is experiment from m%d on %s valid?',id,tmpDate);
            includeSession = input(tmpPrompt);
            if includeSession
                insert(v1inf.Experiment,{id,thisDate,thisDir,thisStimProc}),
            end
        end
    end
end
    
%%
% After adding new session to database from sever, run following code
populate(v1inf.ExpSync),
populate(v1inf.MeanRefImg),
populate(v1inf.Neuron),
populate(v1inf.Target),
populate(v1inf.RandomGratingsExp),
% Manually export randomGratings data for fitting on orchestra
populate(v1inf.StimGratingsData),
populate(v1inf.FiltOverlap),
populate(v1inf.Influence),
populate(v1inf.SelfStim),

%%
% After fitting randomGratings on Orchestra, execute the following code
populate(v1inf.RandomGratingsGP),
populate(v1inf.TuningProps),
populate(v1inf.PairCorrs),
populate(v1inf.TuningDOM),