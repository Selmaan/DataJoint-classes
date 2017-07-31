%{
# Fitted Gaussian Process Models to RandomGratings Experiment
-> v1inf.RandomGratingsExp
-----
gp_fitpar: longblob         # struct with parameters and results of fitting gp model
gp_hyppar: longblob         # cell containing optimized hyperparameters for each neuron's model
gp_ph_hyppar: longblob      # cell containing optimized contrast-phasehyperparameters for each neuron's model
gp_tuning: longblob         # struct with tuning curves extracted from each neuron's model
%}

classdef RandomGratingsGP < dj.Computed

	methods(Access=protected)

		function makeTuples(self, key)
            baseDir = 'F:\DataJoint-FileOut';
            thisFn = sprintf('orchestraStimFitsGP_m%d_%s_jobIndex NaN',...
                key.mouse_id, key.exp_date);
            load(fullfile(baseDir,thisFn)),
            
            gPar = replaceFunctionHandlesWithStrings(gPar);
            key.gp_fitpar = gPar;
            key.gp_hyppar = oH;
            key.gp_ph_hyppar = oPh;
            key.gp_tuning = tR;
            
            self.insert(key)
		end
	end

end

function gPar = replaceFunctionHandlesWithStrings(gPar)

gPar.infFunc = func2str(gPar.infFunc);
gPar.likFunc{1} = func2str(gPar.likFunc{1});
gPar.covFunc{1} = func2str(gPar.covFunc{1});
gPar.covFunc{2}{1} = func2str(gPar.covFunc{2}{1});
for i=1:4
    gPar.covFunc{2}{2}{i}{1} = func2str(gPar.covFunc{2}{2}{i}{1});
    gPar.covFunc{2}{2}{i}{2}{2}{1} = func2str(gPar.covFunc{2}{2}{i}{2}{2}{1});
end

gPar.covFuncPh = func2str(gPar.covFuncPh);
gPar.likFuncPh = func2str(gPar.likFuncPh);

end