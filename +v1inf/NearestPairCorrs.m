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
classdef NearestPairCorrs < dj.Computed
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
            theseTargets = v1inf.Target & key;
            nTargets = length(fetchn(theseTargets,'targ_id'));
            theseTunings = v1inf.TuningProps & key;
            nNeurons = length(fetchn(theseTunings,'neur_id'));
            
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
            
            targetTuning = (proj(theseTargets)*v1inf.TargetedNeurons) * theseTunings;
            [targID,targ_testcorr_tmp,targ_traincorr_tmp,targ_predresp_tmp,targ_residresp_tmp,...
                targ_dirmean_tmp,targ_sfmean_tmp,targ_tfmean_tmp,targ_spdmean_tmp] = fetchn(targetTuning,...
                'targ_id','test_corr','train_corr','pred_resp','resid_resp',...
                'dir_mean','sf_mean','tf_mean','spd_mean','ORDER BY targ_id');
            
            % Handle targets with invalid models by addind to invalid targ
            noFitTargets = find(~isfinite(targ_testcorr_tmp));
            if ~isempty(noFitTargets)
                fprintf('Found %d Invalid Target Model(s) \n',length(noFitTargets)),
                targID = setdiff(targID,targID(noFitTargets));
            end
            invalidTargets = find(~ismember(1:nTargets,targID));
                
            targ_testcorr(targID) = targ_testcorr_tmp(isfinite(targ_testcorr_tmp));
            targ_traincorr(targID) = targ_traincorr_tmp(isfinite(targ_testcorr_tmp));
            targ_predresp(:,targID) = cat(2,targ_predresp_tmp{:});
            targ_residresp(:,targID) = cat(2,targ_residresp_tmp{:});
            targ_dirmean(:,targID) = cat(2,targ_dirmean_tmp{:});
            targ_sfmean(:,targID) = cat(2,targ_sfmean_tmp{:});
            targ_tfmean(:,targID) = cat(2,targ_tfmean_tmp{:});
            targ_spdmean(:,targID) = cat(2,targ_spdmean_tmp{:});
            
            targ_testcorr(invalidTargets) = nan;
            targ_traincorr(invalidTargets) = nan;
            targ_predresp(:,invalidTargets) = nan;
            targ_residresp(:,invalidTargets) = nan;
            targ_dirmean(:,invalidTargets) = nan;
            targ_sfmean(:,invalidTargets) = nan;
            targ_tfmean(:,invalidTargets) = nan;
            targ_spdmean(:,invalidTargets) = nan;
            
            nt_sigcorr = corr(neur_predresp,targ_predresp);
            nt_noicorr = corr(neur_residresp,targ_residresp);
            nt_respcorr = corr(neur_predresp+neur_residresp,targ_predresp+targ_residresp);
            
            neur_dirtune = neur_dirmean(1:180,:)-neur_dirmean(181:360,:);
            neur_oritune = neur_dirmean(1:180,:)+neur_dirmean(181:360,:);
            targ_dirtune = targ_dirmean(1:180,:)-targ_dirmean(181:360,:);
            targ_oritune = targ_dirmean(1:180,:)+targ_dirmean(181:360,:);
            
            nt_dircorr = corr(neur_dirtune,targ_dirtune);
            nt_oricorr = corr(neur_oritune,targ_oritune);
            nt_sfcorr = corr(neur_sfmean,targ_sfmean);
            nt_tfcorr = corr(neur_tfmean,targ_tfmean);
            nt_spdcorr = corr(neur_spdmean,targ_spdmean);
            nt_tuningcorr = corr([neur_dirmean;neur_sfmean;neur_tfmean;neur_spdmean],...
                [targ_dirmean;targ_sfmean;targ_tfmean;targ_spdmean]);
            
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