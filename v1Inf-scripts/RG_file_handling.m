%% Create allResp .mat files for fitting on orchestra

allExps = fetch(v1inf.RandomGratingsExp & (v1inf.ExpType & 'exp_type="multi-contrast"'));
for i=1:length(allExps)
    thisExp = v1inf.RandomGratingsExp & allExps(i);
    allResp = fetch1(thisExp,'rg_resps');
    
    thisFn = sprintf('orchestraStimFitsGP_m%d_%s',...
        allExps(i).mouse_id, allExps(i).exp_date);
    thisDir = fullfile('F:\DataJoint-FileOut',thisFn);
    save(thisDir,'allResp'),
end