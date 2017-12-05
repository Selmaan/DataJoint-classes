%{
# Pre-Computed Shuffled Correlations for each target-neuron pair
-> v1inf.Neuron
-> v1inf.Target
-----
shuf_sigcorr: blob  # Signal Correlation
shuf_noicorr: blob  # Noise Correlation
shuf_respcorr: blob  # Response Correlation
shuf_dircorr: blob  # Direction Correlation
shuf_oricorr: blob  # Orientation Correlation
shuf_sfcorr: blob  # Spatial Frequency Correlation
shuf_tfcorr: blob  # Temporal Frequency Correlation
shuf_spdcorr: blob  # Running Speed Correlation
shuf_tuningcorr: blob  # Tuning Curve Correlation
shuf_neur_testcorr: blob # Neuron Tuning Model Test Correlation
shuf_neur_traincorr: blob # Neuron Tuning Model Train Correlation
shuf_targ_testcorr: blob # Target Tuning Model Test Correlation
shuf_targ_traincorr: blob # Target Tuning Model Train Correlation
%}


%%
classdef ShufflePairCorr < dj.Computed
    properties
        popRel = v1inf.ExpSync;
    end
    
    methods(Access=protected)

		function makeTuples(self, key)
            nShuffles = 1e3;
            
            % Set up Keys and Variables
            theseTargets = v1inf.Target & key;
            nTargets = length(fetchn(theseTargets,'targ_id'));
            theseTunings = v1inf.TuningProps & key;
            nNeurons = length(fetchn(theseTunings,'neur_id'));
            
            key.shuf_sigcorr = nan(nShuffles,1);
            key.shuf_noicorr = nan(nShuffles,1);
            key.shuf_respcorr = nan(nShuffles,1);
            key.shuf_dircorr = nan(nShuffles,1);
            key.shuf_oricorr = nan(nShuffles,1);
            key.shuf_sfcorr = nan(nShuffles,1);
            key.shuf_tfcorr = nan(nShuffles,1);
            key.shuf_spdcorr = nan(nShuffles,1);
            key.shuf_tuningcorr = nan(nShuffles,1);
            key.shuf_neur_testcorr = nan(nShuffles,1);
            key.shuf_neur_traincorr = nan(nShuffles,1);
            key.shuf_targ_testcorr = nan(nShuffles,1);
            key.shuf_targ_traincorr = nan(nShuffles,1);
            
            keys = repmat(key,[nNeurons nTargets]);
            for nNeur = 1:nNeurons
                for nTarg = 1:nTargets
                    keys(nNeur,nTarg).neur_id = nNeur;
                    keys(nNeur,nTarg).targ_id = nTarg;
                end
            end
            
            % Load Neuron Data
            [neurID,neur_testcorr_tmp,neur_traincorr_tmp,neur_predresp_tmp,neur_residresp_tmp,...
                neur_dirmean_tmp,neur_sfmean_tmp,neur_tfmean_tmp,neur_spdmean_tmp] = fetchn(theseTunings,...
                'neur_id','test_corr','train_corr','pred_resp','resid_resp',...
                'dir_mean','sf_mean','tf_mean','spd_mean','ORDER BY neur_id');
            % Handle Neurons with invalid models by addind to invalid neur
            noFitNeurons = find(~isfinite(neur_testcorr_tmp));
            if ~isempty(noFitNeurons)
                fprintf('Found %d Invalid Neuron Model(s) \n',length(noFitNeurons)),
                neurID = setdiff(neurID,noFitNeurons);
            end
            invalidNeurons = find(~ismember(1:nNeurons,neurID));
            
            % Load Target Data
            targetTuning = theseTargets * theseTunings;
            [targID,targ_testcorr_tmp,targ_traincorr_tmp,targ_predresp_tmp,targ_residresp_tmp,...
                targ_dirmean_tmp,targ_sfmean_tmp,targ_tfmean_tmp,targ_spdmean_tmp] = fetchn(targetTuning,...
                'targ_id','test_corr','train_corr','pred_resp','resid_resp',...
                'dir_mean','sf_mean','tf_mean','spd_mean','ORDER BY targ_id');            
            % Handle targets with invalid models by addind to invalid targ
            noFitTargets = find(~isfinite(targ_testcorr_tmp));
            if ~isempty(noFitTargets)
                fprintf('Found %d Invalid Target Model(s) \n',length(noFitTargets)),
                targID = setdiff(targID,noFitTargets);
            end
            invalidTargets = find(~ismember(1:nTargets,targID));
            
            % Perform Shuffles
            for iShuffle = 1:nShuffles
                if mod(iShuffle,100)==1
                    iShuffle,
                end
                % Permute Neuron IDs
                neurID = neurID(randperm(length(neurID)));
                % Permute Target IDs
                targID = targID(randperm(length(targID)));
                % Assemble neuron tunings
                neur_testcorr(neurID) = neur_testcorr_tmp(isfinite(neur_testcorr_tmp));
                neur_traincorr(neurID) = neur_traincorr_tmp(isfinite(neur_testcorr_tmp));
                neur_predresp(:,neurID) = cat(2,neur_predresp_tmp{:});
                neur_residresp(:,neurID) = cat(2,neur_residresp_tmp{:});
                neur_dirmean(:,neurID) = cat(2,neur_dirmean_tmp{:});
                neur_sfmean(:,neurID) = cat(2,neur_sfmean_tmp{:});
                neur_tfmean(:,neurID) = cat(2,neur_tfmean_tmp{:});
                neur_spdmean(:,neurID) = cat(2,neur_spdmean_tmp{:});
                % Set invalid neuron tunings to NAN
                neur_testcorr(invalidNeurons) = nan;
                neur_traincorr(invalidNeurons) = nan;
                neur_predresp(:,invalidNeurons) = nan;
                neur_residresp(:,invalidNeurons) = nan;
                neur_dirmean(:,invalidNeurons) = nan;
                neur_sfmean(:,invalidNeurons) = nan;
                neur_tfmean(:,invalidNeurons) = nan;
                neur_spdmean(:,invalidNeurons) = nan;
                
                % Assemble Target Tunings
                targ_testcorr(targID) = targ_testcorr_tmp(isfinite(targ_testcorr_tmp));
                targ_traincorr(targID) = targ_traincorr_tmp(isfinite(targ_testcorr_tmp));
                targ_predresp(:,targID) = cat(2,targ_predresp_tmp{:});
                targ_residresp(:,targID) = cat(2,targ_residresp_tmp{:});
                targ_dirmean(:,targID) = cat(2,targ_dirmean_tmp{:});
                targ_sfmean(:,targID) = cat(2,targ_sfmean_tmp{:});
                targ_tfmean(:,targID) = cat(2,targ_tfmean_tmp{:});
                targ_spdmean(:,targID) = cat(2,targ_spdmean_tmp{:});
                % Set invalid Targets to NAN
                targ_testcorr(invalidTargets) = nan;
                targ_traincorr(invalidTargets) = nan;
                targ_predresp(:,invalidTargets) = nan;
                targ_residresp(:,invalidTargets) = nan;
                targ_dirmean(:,invalidTargets) = nan;
                targ_sfmean(:,invalidTargets) = nan;
                targ_tfmean(:,invalidTargets) = nan;
                targ_spdmean(:,invalidTargets) = nan;
            
                % Compute tuning correlations
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
            
                % Create Shuffle Key
                for nNeur = 1:size(keys,1)
                    for nTarg = 1:size(keys,2)
                        keys(nNeur,nTarg).shuf_sigcorr(iShuffle) = nt_sigcorr(nNeur,nTarg);
                        keys(nNeur,nTarg).shuf_noicorr(iShuffle) = nt_noicorr(nNeur,nTarg);
                        keys(nNeur,nTarg).shuf_respcorr(iShuffle) = nt_respcorr(nNeur,nTarg);
                        keys(nNeur,nTarg).shuf_dircorr(iShuffle) = nt_dircorr(nNeur,nTarg);
                        keys(nNeur,nTarg).shuf_oricorr(iShuffle) = nt_oricorr(nNeur,nTarg);
                        keys(nNeur,nTarg).shuf_sfcorr(iShuffle) = nt_sfcorr(nNeur,nTarg);
                        keys(nNeur,nTarg).shuf_tfcorr(iShuffle) = nt_tfcorr(nNeur,nTarg);
                        keys(nNeur,nTarg).shuf_spdcorr(iShuffle) = nt_spdcorr(nNeur,nTarg);
                        keys(nNeur,nTarg).shuf_tuningcorr(iShuffle) = nt_tuningcorr(nNeur,nTarg);
                        keys(nNeur,nTarg).shuf_neur_testcorr(iShuffle) = neur_testcorr(nNeur);
                        keys(nNeur,nTarg).shuf_neur_traincorr(iShuffle) = neur_traincorr(nNeur);
                        keys(nNeur,nTarg).shuf_targ_testcorr(iShuffle) = targ_testcorr(nTarg);
                        keys(nNeur,nTarg).shuf_targ_traincorr(iShuffle) = targ_traincorr(nTarg);
                    end
                end
            end
            
            if numel(keys)<1e4
                insert(self,keys(:)),
            else
                insert(self,keys(1:1e4)),
                insert(self,keys(1e4+1:end)),
            end
        end
    end
end