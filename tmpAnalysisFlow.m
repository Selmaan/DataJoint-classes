%%
validStim = v1inf.SelfStim & 'self_stim>.05';
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
[iD,iC,rM,rV] = fetchn(((v1inf.Influence & validFilt) & validStim) - v1inf.Target,...
    'inf_dist','inf_naivecorr','inf_rawmu','inf_rawvar');
iC(isnan(iC)) = 0;

%%
tmpInd = iD>=25;

tmpDist = iD(tmpInd);
tmpInf = rM(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));

dBin = nan(650,1);
for i=1:650
    dBin(i) = median(tmpInf(tmpDist>i & tmpDist<i+40));
end

hold on,
plot(21:670,dBin,'linewidth',2)
axis tight,
xlabel('Neuron-Target Distance (um)')
%%
tmpInd = iD>70;
tmpCorr = iC(tmpInd);
% tmpInf = rM(tmpInd);
tmpInf = rM(tmpInd)./sqrt(rV(tmpInd));
% cVals = [-inf,linspace(min(tmpCorr)-.01,max(tmpCorr)+.01,13),inf];
cVals = linspace(-.4,.7,20);
% cVals = [-1, linspace(-.2,.4,7), 1];

cBin = nan(length(cVals)-1,1);
fracInBin = nan(length(cVals)-1,1);
for i=1:length(cVals)-1
    cBin(i) = mean(tmpInf(tmpCorr>=cVals(i) & tmpCorr<cVals(i+1)));
    fracInBin(i) = mean((tmpCorr>=cVals(i) & tmpCorr<cVals(i+1)));
end

figure,plot(cVals(1:end-1) + mode(diff(cVals)),cBin,'linewidth',2)
figure,semilogy(cVals(1:end-1) + mode(diff(cVals)),fracInBin,'linewidth',2)
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