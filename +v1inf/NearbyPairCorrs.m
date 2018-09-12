%{
# Pre-Computed Correlations for each target-neuron pair (using nearest neuron to target instead of manual match)
-> v1inf.Neuron
-> v1inf.Target
-----
inf_sigcorr=NULL: double  # Signal Correlation
inf_noicorr=NULL: double  # Noise Correlation
inf_respcorr=NULL: double  # Response Correlation
inf_dircorr=NULL: double  # Direction Correlation
inf_oricorr=NULL: double  # Orientation Correlation
inf_sfcorr=NULL: double  # Spatial Frequency Correlation
inf_tfcorr=NULL: double  # Temporal Frequency Correlation
inf_spdcorr=NULL: double  # Running Speed Correlation
inf_tuningcorr=NULL: double  # Tuning Curve Correlation
inf_neur_testcorr=NULL: double # Neuron Tuning Model Test Correlation
inf_neur_traincorr=NULL: double # Neuron Tuning Model Train Correlation
inf_targ_testcorr=NULL: double # Target Tuning Model Test Correlation
inf_targ_traincorr=NULL: double # Target Tuning Model Train Correlation
%}


%%
classdef NearbyPairCorrs < dj.Computed
    properties
        popRel = v1inf.RandomGratingsGP;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            thisExpt = fetch1(v1inf.ExpType & key,'exp_type');
            if contains(thisExpt,'Monitor-Off')
                fprintf('Skipping Random Gratings for Monitor Off Experiment \n');
                return
            end
            theseTargets = fetch(v1inf.NearbyNeurons & key,'*');
            nTargets = length(theseTargets);
            theseTunings = v1inf.TuningProps & key;
            nNeurons = length(fetchn(theseTunings,'neur_id'));
            
            nt_sigcorr = nan(nNeurons,nTargets);
            nt_noicorr = nan(nNeurons,nTargets);
            nt_respcorr = nan(nNeurons,nTargets);
            nt_dircorr = nan(nNeurons,nTargets);
            nt_oricorr = nan(nNeurons,nTargets);
            nt_sfcorr = nan(nNeurons,nTargets);
            nt_tfcorr = nan(nNeurons,nTargets);
            nt_spdcorr = nan(nNeurons,nTargets);
            nt_tuningcorr = nan(nNeurons,nTargets);
            
            [neurID,neur_testcorr_tmp,neur_traincorr_tmp,neur_predresp_tmp,neur_residresp_tmp,...
                neur_dirmean_tmp,neur_sfmean_tmp,neur_tfmean_tmp,neur_spdmean_tmp] = fetchn(theseTunings,...
                'neur_id','test_corr','train_corr','pred_resp','resid_resp',...
                'dir_mean','sf_mean','tf_mean','spd_mean','ORDER BY neur_id');
            % Handle Neurons with invalid models by addind to invalid neur
            noFitNeurons = find(~isfinite(neur_testcorr_tmp));
            if ~isempty(noFitNeurons)
                fprintf('Found %d Invalid Neuron Model(s) \n',length(noFitNeurons)),
                neurID = setdiff(neurID,neurID(noFitNeurons));
            end
            invalidNeurons = find(~ismember(1:nNeurons,neurID));
            
            neur_testcorr(neurID) = neur_testcorr_tmp(isfinite(neur_testcorr_tmp));
            neur_traincorr(neurID) = neur_traincorr_tmp(isfinite(neur_testcorr_tmp));
            neur_predresp(:,neurID) = cat(2,neur_predresp_tmp{:});
            neur_residresp(:,neurID) = cat(2,neur_residresp_tmp{:});
            neur_dirmean(:,neurID) = cat(2,neur_dirmean_tmp{:});
            neur_sfmean(:,neurID) = cat(2,neur_sfmean_tmp{:});
            neur_tfmean(:,neurID) = cat(2,neur_tfmean_tmp{:});
            neur_spdmean(:,neurID) = cat(2,neur_spdmean_tmp{:});
            
            neur_testcorr(invalidNeurons) = nan;
            neur_traincorr(invalidNeurons) = nan;
            neur_predresp(:,invalidNeurons) = nan;
            neur_residresp(:,invalidNeurons) = nan;
            neur_dirmean(:,invalidNeurons) = nan;
            neur_sfmean(:,invalidNeurons) = nan;
            neur_tfmean(:,invalidNeurons) = nan;
            neur_spdmean(:,invalidNeurons) = nan;
            
            neur_dirtune = neur_dirmean(1:180,:)-neur_dirmean(181:360,:);
            neur_oritune = neur_dirmean(1:180,:)+neur_dirmean(181:360,:);
            
            % Loop over each target site individually, and grab all nearby
            % neurons for comparison with other neurons (excludes <25um)
            for n = 1:nTargets
                nNearbyNeur = length(theseTargets(n).neur_ids);
                neurKeys = repmat(rmfield(theseTargets(n),'neur_ids'),nNearbyNeur,1);
                for i = 1:nNearbyNeur
                    neurKeys(i).neur_id = theseTargets(n).neur_ids(i);
                end
                % Do not compare a neuron to itself!
                validComparisons = setdiff(1:nNeurons, theseTargets(n).neur_ids);
                
                targetTuning = v1inf.TuningProps & neurKeys;
                [targ_testcorr_tmp,targ_traincorr_tmp,targ_predresp_tmp,targ_residresp_tmp,...
                    targ_dirmean_tmp,targ_sfmean_tmp,targ_tfmean_tmp,targ_spdmean_tmp] = fetchn(targetTuning,...
                    'test_corr','train_corr','pred_resp','resid_resp',...
                    'dir_mean','sf_mean','tf_mean','spd_mean');
                targ_predresp = cell2mat(targ_predresp_tmp');
                targ_residresp = cell2mat(targ_residresp_tmp');
                targ_dirmean = cell2mat(targ_dirmean_tmp');
                targ_sfmean = cell2mat(targ_sfmean_tmp');
                targ_tfmean = cell2mat(targ_tfmean_tmp');
                targ_spdmean = cell2mat(targ_spdmean_tmp');
                targ_dirtune = targ_dirmean(1:180,:)-targ_dirmean(181:360,:);
                targ_oritune = targ_dirmean(1:180,:)+targ_dirmean(181:360,:);
                
                targ_testcorr(n) = nanmean(targ_testcorr_tmp);
                targ_traincorr(n) = nanmean(targ_traincorr_tmp);
                nt_sigcorr(validComparisons,n) = nanmean(corr(neur_predresp(:,validComparisons),targ_predresp),2);
                nt_noicorr(validComparisons,n) = nanmean(corr(neur_residresp(:,validComparisons),targ_residresp),2);
                nt_respcorr(validComparisons,n) = nanmean(corr(neur_predresp(:,validComparisons)+neur_residresp(:,validComparisons),...
                    targ_predresp+targ_residresp),2);
                
                nt_dircorr(validComparisons,n) = nanmean(corr(neur_dirtune(:,validComparisons),targ_dirtune),2);
                nt_oricorr(validComparisons,n) = nanmean(corr(neur_oritune(:,validComparisons),targ_oritune),2);
                nt_sfcorr(validComparisons,n) = nanmean(corr(neur_sfmean(:,validComparisons),targ_sfmean),2);
                nt_tfcorr(validComparisons,n) = nanmean(corr(neur_tfmean(:,validComparisons),targ_tfmean),2);
                nt_spdcorr(validComparisons,n) = nanmean(corr(neur_spdmean(:,validComparisons),targ_spdmean),2);
                nt_tuningcorr(validComparisons,n) = nanmean(corr([neur_dirmean(:,validComparisons);neur_sfmean(:,validComparisons);neur_tfmean(:,validComparisons);neur_spdmean(:,validComparisons)],...
                    [targ_dirmean;targ_sfmean;targ_tfmean;targ_spdmean]),2);
            
            end
            
            keys = repmat(key,size(nt_sigcorr));
            for nNeur = 1:size(keys,1)
                for nTarg = 1:size(keys,2)
                    keys(nNeur,nTarg).neur_id = nNeur;
                    keys(nNeur,nTarg).targ_id = nTarg;
                    keys(nNeur,nTarg).inf_sigcorr = nt_sigcorr(nNeur,nTarg);
                    keys(nNeur,nTarg).inf_noicorr = nt_noicorr(nNeur,nTarg);
                    keys(nNeur,nTarg).inf_respcorr = nt_respcorr(nNeur,nTarg);
                    keys(nNeur,nTarg).inf_dircorr = nt_dircorr(nNeur,nTarg);
                    keys(nNeur,nTarg).inf_oricorr = nt_oricorr(nNeur,nTarg);
                    keys(nNeur,nTarg).inf_sfcorr = nt_sfcorr(nNeur,nTarg);
                    keys(nNeur,nTarg).inf_tfcorr = nt_tfcorr(nNeur,nTarg);
                    keys(nNeur,nTarg).inf_spdcorr = nt_spdcorr(nNeur,nTarg);
                    keys(nNeur,nTarg).inf_tuningcorr = nt_tuningcorr(nNeur,nTarg);
                    keys(nNeur,nTarg).inf_neur_testcorr = neur_testcorr(nNeur);
                    keys(nNeur,nTarg).inf_neur_traincorr = neur_traincorr(nNeur);
                    keys(nNeur,nTarg).inf_targ_testcorr = targ_testcorr(nTarg);
                    keys(nNeur,nTarg).inf_targ_traincorr = targ_traincorr(nTarg);
                end
            end
                    
            insert(self,keys),
            
        end
    end
end