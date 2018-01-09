function [yVal,xVal,yStd] = prctileSmoothingBootstrap(X,Y,halfWidth,nBootstraps)

prctileEdges = prctile(X,0:100);

xVal = nan(101,1);
yVal = nan(101,1);
yStd = nan(101,1);
for prctileBin = 1:101
    minVal = prctileEdges(max(prctileBin-halfWidth,1));
    maxVal = prctileEdges(min(prctileBin+halfWidth,101));
    binInd = X<=maxVal & X>=minVal;
    xVal(prctileBin) = mean(X(binInd));
    theseVals = Y(binInd);
    yVal(prctileBin) = mean(theseVals);
    nBinCounts = length(theseVals);
    thisBootstrap = nan(nBootstraps,1);
    parfor nBootstrap = 1:nBootstraps
        theseInd = randi(nBinCounts,nBinCounts,1);
        thisBootstrap(nBootstrap) = mean(theseVals(theseInd));
    end
    yStd(prctileBin) = std(thisBootstrap);
end