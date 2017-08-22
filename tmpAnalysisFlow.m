%%
validStim = v1inf.SelfStim & 'self_stim>5';
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
[iD,iC,rM,rV,sN,sP,keys] = fetchn(((v1inf.Influence & validFilt) & validStim) - v1inf.Target,...
    'inf_dist','inf_naivecorr','inf_regmu','inf_regvar','inf_shuf_n','inf_shuf_p');
iC(isnan(iC)) = 0;

%%
tmpInd = iD>=25;

tmpDist = iD(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
tmpInf = sN(tmpInd);

dBin = nan(650,1);
for i=1:650
    dBin(i) = mean(tmpInf(tmpDist>i & tmpDist<i+50));
end

hold on,
plot(26:675,dBin,'linewidth',2)
axis tight,
xlabel('Neuron-Target Distance (um)')
%%
tmpInd = iD>30;
tmpCorr = iC(tmpInd);
tmpInf = sN(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
% cVals = [-inf,linspace(min(tmpCorr)-.01,max(tmpCorr)+.01,13),inf];
% cVals = linspace(-.4,.7,10);
cVals = [-1, linspace(-.15,.3,7), 1];

cBin = nan(length(cVals)-1,1);
xBin = nan(length(cVals)-1,1);
numInBin = nan(length(cVals)-1,1);
for i=1:length(cVals)-1
    validInds=tmpCorr>=cVals(i) & tmpCorr<cVals(i+1);
    cBin(i) = nanmean(tmpInf(validInds));
    xBin(i) = nanmean(tmpCorr(validInds));
    numInBin(i) = sum((tmpCorr>=cVals(i) & tmpCorr<cVals(i+1)));
end

hold on,plot(xBin,cBin,'linewidth',2)
hold on,plot(xBin,cBin,'*','markersize',10),
% figure,semilogy(xBin,numInBin,'linewidth',2)
% hold on,plot(xBin,numInBin,'*','markersize',10),
%%
% distBins = [0,30,70,110,170,250,400,550,inf];
distBins = [0,50,250,400,inf];
corrBins = [-inf, -.15, -.05, 0.01, .05, .2, inf];

infBins = [];
numInBin = [];
tmpInf = rM./sqrt(rV);
for distBin = 1:length(distBins)-1
    distInd = iD>distBins(distBin) & iD < distBins(distBin+1);
    for corrBin = 1:length(corrBins)-1
        corrInd = iC>corrBins(corrBin) & iC < corrBins(corrBin+1);
        infBins(corrBin,distBin) = mean(tmpInf(distInd & corrInd));
        numInBin(corrBin,distBin) = sum(distInd & corrInd);
    end
end