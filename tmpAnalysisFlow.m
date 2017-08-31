%%
validBin = v1inf.InfluenceBin & 'bin_id=12';
validStim = v1inf.SelfStim & 'self_stim>5';
validFilt = v1inf.FiltOverlap - 'filt_overlap>0';
[iD,iC,rM,rV,sN,sP] = fetchn(((validBin & validFilt) & validStim) - v1inf.Target,...
    'inf_dist','inf_naivecorr','inf_regmu','inf_regvar','inf_shuf_n','inf_shuf_p');
iC(isnan(iC)) = 0;

%%
tmpInd = iD>=25;

tmpDist = iD(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
tmpInf = sP(tmpInd);

dBin = nan(650,1);
for i=1:650
    dBin(i) = nanmean(tmpInf(tmpDist>i & tmpDist<i+50));
end

hold on,
plot(26:675,dBin,'linewidth',2)
axis tight,
xlabel('Neuron-Target Distance (um)')
%%
tmpInd = iD>30;
tmpCorr = iC(tmpInd);
tmpInf = sP(tmpInd);
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

%% Distance Only GP Model

hyp = struct();
likFunc = {@likGauss};
hyp.lik = 0;
% covFunc = {@covSEiso};
% hyp.cov = [3;-1];
covFunc = {@covRQiso};
hyp.cov = [3;-1;0];
meanFunc = {@meanConst};
hyp.mean = 0;
infFunc = @infGaussLik;

xu = linspace(25,700,30)'; covFunc = {'apxSparse',covFunc,xu};      % inducing points
infApx  = @(varargin) infFunc(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
hyp.xu = xu;

tmpInd = iD>=30;
tmpDist = iD(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
tmpInf = sN(tmpInd);
optHyp = minimize(hyp,@gp,-75,infApx,meanFunc,covFunc,likFunc,tmpDist,tmpInf);

xDisp = (25:650)';
[~,~,yMu,yVar] = gp(optHyp,infApx,meanFunc,covFunc,likFunc,tmpDist,tmpInf,xDisp);
figure,plot(xDisp,yMu,'k')
hold on,plot(xDisp,yMu+sqrt(yVar),'r')
hold on,plot(xDisp,yMu-sqrt(yVar),'r')

%% Correlation Only GP Model
tmpInd = iD>=30;
tmpCorr = iC(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
% tmpCorr(isnan(tmpInf))=[]; tmpInf(isnan(tmpInf))=[];
tmpInf = sN(tmpInd);

hyp = struct();
likFunc = {@likGauss};
hyp.lik = 0;
% covFunc = {@covSEiso};
% hyp.cov = [-1.5;-1];
covFunc = {@covRQiso};
hyp.cov = [-2;-1;0];
meanFunc = {@meanConst};
hyp.mean = 0;
infFunc = @infGaussLik;

xu = linspace(min(tmpCorr),max(tmpCorr),20)'; covFunc = {'apxSparse',covFunc,xu};      % inducing points
infApx  = @(varargin) infFunc(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
hyp.xu = xu;

optHyp = minimize(hyp,@gp,-75,infApx,meanFunc,covFunc,likFunc,tmpCorr,tmpInf);

xDisp = linspace(min(tmpCorr),max(tmpCorr),100)';
[~,~,yMu,yVar] = gp(optHyp,infApx,meanFunc,covFunc,likFunc,tmpCorr,tmpInf,xDisp);
figure,plot(xDisp,yMu,'k')
hold on,plot(xDisp,yMu+sqrt(yVar),'r')
hold on,plot(xDisp,yMu-sqrt(yVar),'r')

%% 2-d regression, additive covariance
tmpInd = iD>=30;
tmpDist = iD(tmpInd);
tmpCorr = iC(tmpInd);
% tmpInf = (rM(tmpInd)./sqrt(rV(tmpInd)));
% tmpCorr(isnan(tmpInf))=[]; tmpDist(isnan(tmpInf)) = []; tmpInf(isnan(tmpInf))=[];
tmpInf = sN(tmpInd);

hyp = struct();
likFunc = {@likGauss};
hyp.lik = 0;
meanFunc = {@meanConst};
hyp.mean = 0;
infFunc = @infGaussLik;
% covFun1 = {@covMask, {[1 0], {@covRQ,'iso',[]}}}';
% covFun2 = {@covMask, {[0 1], {@covRQ,'iso',[]}}}';
% covFunc = {@covScale {@covSum, {covFun1, covFun2}}};
% hyp.cov = [3;0;-2;0;0];
covFunc = {@covADD,{[1,2],'covRQiso'}};
hyp.cov = [3;0;0;-2;0;0;0;0];

xu_d = linspace(25,700,20)'; 
xu_c = linspace(min(tmpCorr),max(tmpCorr),20)';
[xu_d,xu_c] = ndgrid(xu_d,xu_c);
xu = [xu_d(:), xu_c(:)];

covFunc = {'apxSparse',covFunc,xu};      % inducing points
infApx  = @(varargin) infFunc(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
hyp.xu = xu;

optHyp = minimize(hyp,@gp,-100,...
    infApx,meanFunc,covFunc,likFunc,[tmpDist,tmpCorr],tmpInf);

[disp_d, disp_c] = ndgrid(linspace(25,650,100)',...
    linspace(min(tmpCorr),max(tmpCorr),99)'); 
xDisp = [disp_d(:), disp_c(:)];
[~,~,yMu,yVar] = gp(optHyp,infApx,meanFunc,covFunc,likFunc,[tmpDist,tmpCorr],tmpInf,xDisp);
yMu = reshape(yMu,100,99); yVar = reshape(yVar,100,99);
figure,imagesc(disp_c(1,:),disp_d(:,1),yMu)
colorbar
xlabel('Neuron-Target Trace Correlation')
ylabel('Neuron Target Distance (um)')

figure,plot(disp_c(1,:),yMu([4,18,50,75,95],:)','linewidth',2)
axis tight
xlabel('Neuron-Target Correlation')
ylabel('Mean Influence')

figure,plot(disp_d(:,1),yMu(:,[5,25,50,75,95]),'linewidth',2)
legend('-.3','-.1','.06','.32','.58')
axis tight
xlabel('Neuron-Target Distance (um)')
ylabel('Mean Influence')
title('Additive Covariance'),

