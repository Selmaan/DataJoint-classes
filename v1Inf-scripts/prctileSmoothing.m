function [yVal,xVal] = prctileSmoothing(X,Y,halfWidth)

prctileEdges = prctile(X,0:100);

xVal = nan(101,1);
yVal = nan(101,1);
for prctileBin = 1:101
    minVal = prctileEdges(max(prctileBin-halfWidth,1));
    maxVal = prctileEdges(min(prctileBin+halfWidth,101));
    binInd = X<=maxVal & X>=minVal;
    xVal(prctileBin) = mean(X(binInd));
    yVal(prctileBin) = mean(Y(binInd));
end