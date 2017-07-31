%% Create allResp .mat files for fitting on orchestra

allExps = fetch(v1inf.ExpSync)
for i=1:length(allExps)
    q1 = sprintf('mouse_id = %d', allExps(i).mouse_id);
    q2 = sprintf('exp_date = "%s"',allExps(i).exp_date);
    thisQ = v1inf.RandomGratingsExp & q1 & q2,
    allResp = fetch1(thisQ,'rg_resps');
    
    thisFn = sprintf('orchestraStimFitsGP_m%d_%s',...
        allExps(i).mouse_id, allExps(i).exp_date);
    thisDir = fullfile('F:\DataJoint-FileOut',thisFn);
    save(thisDir,'allResp'),
end