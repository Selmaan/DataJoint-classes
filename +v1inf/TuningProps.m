%{
# Extracted Tuning Properties for indv neuron from GP model
-> v1inf.RandomGratingsGP
-> v1inf.Neuron
-----
test_corr=NULL : double      # corr btw predicted and actual grating resp on held-out data
train_corr=NULL : double     # corr btw predicted and actual grating resp on training data
pred_resp=NULL : blob        # predicted gratings resp on held-out data
resid_resp=NULL : blob       # actual minus predicted resp on held-out data
phase_resp=NULL : blob       # grating contrast-phase response
lik_hyp=NULL : double        # log(noise std) from fitted grating GP model
dir_hyp=NULL : double        # log(dir length) ...
sf_hyp=NULL : double         # log(SF length) ...
tf_hyp=NULL : double         # log(TF length) ...
spd_hyp=NULL : double        # log(run spd length) ...
cov_hyp=NULL : double        # log(cov scale) ...
dir_mean=NULL : blob         # Direction Tuning, mean of latent function
sf_mean=NULL : blob          # SF Tuning, ...
tf_mean=NULL : blob          # TF Tuning, ...
spd_mean=NULL : blob         # Running spd Tuning, ...
dir_std=NULL : blob          # Direction Tuning, std of latent function
sf_std=NULL : blob           # SF Tuning, ...
tf_std=NULL : blob           # TF Tuning, ...
spd_std=NULL : blob          # Running spd Tuning, ...
%}

classdef TuningProps < dj.Computed
    
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
            theseNeurons = fetch(v1inf.Neuron & key);
            thisGP = fetch(v1inf.RandomGratingsGP & key,'*');
            nNeurons = length(theseNeurons);
            if length(thisGP.gp_hyppar) ~= nNeurons
                error('Inconsistent Data Structure!'),
            end
            
            keys = repmat(key,nNeurons,1);
            for nNeuron = 1:nNeurons
                keys(nNeuron).neur_id = nNeuron;
                if ~isempty(thisGP.gp_hyppar(nNeuron).lik)
                    keys(nNeuron).test_corr = thisGP.gp_fitpar.predCorr(nNeuron,1);
                    keys(nNeuron).train_corr = thisGP.gp_fitpar.predCorr(nNeuron,2);
                    keys(nNeuron).pred_resp = thisGP.gp_fitpar.testPreds(:,nNeuron);
                    keys(nNeuron).resid_resp = thisGP.gp_fitpar.Y(:,nNeuron) - thisGP.gp_fitpar.testPreds(:,nNeuron);
                    keys(nNeuron).phase_resp = thisGP.gp_fitpar.yPh(:,nNeuron);
                    keys(nNeuron).lik_hyp = thisGP.gp_hyppar(nNeuron).lik;
                    keys(nNeuron).dir_hyp = thisGP.gp_hyppar(nNeuron).cov(1);
                    keys(nNeuron).sf_hyp = thisGP.gp_hyppar(nNeuron).cov(2);
                    keys(nNeuron).tf_hyp = thisGP.gp_hyppar(nNeuron).cov(3);
                    keys(nNeuron).spd_hyp = thisGP.gp_hyppar(nNeuron).cov(4);
                    keys(nNeuron).cov_hyp = thisGP.gp_hyppar(nNeuron).cov(5);
                    keys(nNeuron).dir_mean = thisGP.gp_tuning.neurTuning(:,nNeuron,1);
                    keys(nNeuron).sf_mean = thisGP.gp_tuning.neurTuning(:,nNeuron,2);
                    keys(nNeuron).tf_mean = thisGP.gp_tuning.neurTuning(:,nNeuron,3);
                    keys(nNeuron).spd_mean = thisGP.gp_tuning.neurTuning(:,nNeuron,4);
                    keys(nNeuron).dir_std = sqrt(thisGP.gp_tuning.neurTuning(:,nNeuron,5));
                    keys(nNeuron).sf_std = sqrt(thisGP.gp_tuning.neurTuning(:,nNeuron,6));
                    keys(nNeuron).tf_std = sqrt(thisGP.gp_tuning.neurTuning(:,nNeuron,7));
                    keys(nNeuron).spd_std = sqrt(thisGP.gp_tuning.neurTuning(:,nNeuron,8));
                else
                    warning('Missing Neuron %d',nNeuron),
                end
            end
            
            self.insert(keys)
		end
	end

end
