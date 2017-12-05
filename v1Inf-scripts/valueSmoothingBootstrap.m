function [yVal,xVal,yStd,yStd2] = valueSmoothingBootstrap(X,Y,centers,halfWidth,nBootstraps)

nBins = length(centers);
xVal = nan(nBins,1);
yVal = nan(nBins,1);
yStd = nan(nBins,1);
yStd2 = nan(nBins,1);
for nBin = 1:nBins
    minVal = centers(nBin)-halfWidth;
    maxVal = centers(nBin)+halfWidth;
    binInd = X<=maxVal & X>=minVal;
    xVal(nBin) = mean(X(binInd));
    theseVals = Y(binInd);
    yVal(nBin) = mean(theseVals);
    nBinCounts = length(theseVals);
    thisBootstrap = nan(nBootstraps,1);
    parfor nBootstrap = 1:nBootstraps
        theseInd = randi(nBinCounts,nBinCounts,1);
        thisBootstrap(nBootstrap) = mean(theseVals(theseInd));
    end
    yStd(nBin) = std(thisBootstrap);
    yStd2(nBin) = (mad(theseVals)*1.253)/sqrt(nBinCounts);
end