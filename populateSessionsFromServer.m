mouseIDs = fetchn(v1inf.Mouse,'mouse_id');

for nID = 1:length(mouseIDs)
    id = mouseIDs(nID);
    baseDir = sprintf('%s%d%s','Z:\HarveyLab\Selmaan\V1-RF\m',id);
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
    