%{
# Influence Measurement Directed Pair (Regression-based)
-> v1inf.Neuron
-> v1inf.Target
-----
infcont_dist: double            # Pairwise distance (in um)
infcont_rawmu=NULL: double      # Mean influence of stim period response
infcont_rawvar=NULL: double     # Variance of stim period response estimate
%}

classdef InfluenceControl < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            theseNeurons = v1inf.Neuron & key;
            theseTargets = v1inf.Target & key & 'targ_label=0';
            theseStimInfo = v1inf.StimGratingsData & key;
            
            [targ_id, targ_xc, targ_yc] = fetchn(...
                theseTargets,'targ_id','targ_xc','targ_yc','ORDER BY targ_id');
            if isempty(targ_id)
                warning('Experiment on %s contains no control stimulation',key.exp_date),
                return
            end
            
            [neur_xc, neur_yc] = fetchn(...
                theseNeurons,'neur_xc','neur_yc','ORDER BY neur_id');
            stimData = fetch(theseStimInfo,'*');
            
            infDist = sqrt((neur_xc-targ_xc').^2 + (neur_yc-targ_yc').^2);
            [rawMu,rawVar] = computeInfControlRegression(stimData,targ_id);
            
            keys = repmat(key,size(infDist));
            for nNeuron = 1:size(infDist,1)
                for nTarg = 1:size(infDist,2)
                    keys(nNeuron,nTarg).neur_id = nNeuron;
                    keys(nNeuron,nTarg).targ_id = targ_id(nTarg);
                    keys(nNeuron,nTarg).infcont_dist = infDist(nNeuron,nTarg);
                    keys(nNeuron,nTarg).infcont_rawmu = rawMu(nNeuron,nTarg);
                    keys(nNeuron,nTarg).infcont_rawvar = rawVar(nNeuron,nTarg);
                end
            end
            
            insert(self,keys),
            
        end
    end
end


function [rawMu,rawVar] = ...
    computeInfControlRegression(stimData,targ_id)

nStim = length(targ_id);

allDir = unique(stimData.stim_vis_dir);
X = zeros(size(stimData.stim_de_resp,1),nStim+length(allDir));
for iStim = 1:nStim
    X(:,iStim) = (stimData.stim_targ_id == targ_id(iStim));
end

for iDir = 1:length(allDir)
    X(:,nStim+iDir) = (stimData.stim_vis_dir == allDir(iDir));
end

logSpd = log(stimData.stim_mv_spd);
mvBins = 0:25:100;
for mvBin = 1:length(mvBins)-1
    binMin = prctile(logSpd,mvBins(mvBin));
    binMax = prctile(logSpd,mvBins(mvBin+1));
    X(:,nStim+length(allDir)+mvBin) = ...
        (logSpd <= binMax) & (logSpd > binMin);
end
X(:,nStim+length(allDir)+1) = []; % Use lowest movement bin as 'intercept'

validTrials = sum(X(:,1:nStim),2) == 1;
X = X(validTrials,:);
yRaw = stimData.stim_de_resp(validTrials,:);

fitMu_raw = nan(size(X,2)-nStim+1,size(yRaw,2),nStim);
fitVar_raw = nan(size(X,2)-nStim+1,size(yRaw,2),nStim);

for thisControl = 1:nStim
    thisX = X;
    thisX(:,(1:nStim)~=thisControl) = [];
    blmObj = diffuseblm(size(thisX,2),'Intercept',false);
    
    for n=1:size(yRaw,2)
        [~,tmpMu,tmpCov] = estimate(blmObj,thisX,yRaw(:,n),'Display',false);
        tmpVar = diag(tmpCov);
        fitMu_raw(:,n,thisControl) = tmpMu;
        fitVar_raw(:,n,thisControl) = tmpVar;
    end
end

%% Reformat Fitted parameters
rawMu = nan(size(yRaw,2),nStim);
rawVar = nan(size(yRaw,2),nStim);
for thisControl = 1:nStim
    rawMu(:,thisControl) = fitMu_raw(1,:,thisControl);
    rawVar(:,thisControl) = fitVar_raw(1,:,thisControl);
end

end