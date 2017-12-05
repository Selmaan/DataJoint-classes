function gpFit = gp1d_smooth(X,Y)

gpFit = struct();

hyp = struct();
likFunc = {@likGauss};
hyp.lik = log(std(Y));
covFunc = {@covSEiso};
hyp.cov = [log(std(X));0];
meanFunc = {@meanConst};
hyp.mean = 0;
infFunc = @infGaussLik;

xu = linspace(prctile(X,0.1),prctile(X,99.9),50)'; 
covFunc = {'apxSparse',covFunc,xu};      % inducing points
infApx  = @(varargin) infFunc(varargin{:},struct('s',1));           % VFE, opt.s = 0

optHyp = minimize(hyp,@gp,-100,infApx,meanFunc,covFunc,likFunc,X,Y);

xDisp = linspace(prctile(X,0.01),prctile(X,99.99),1e3)';
[~,~,yMu,yVar] = gp(optHyp,infApx,meanFunc,covFunc,likFunc,X,Y,xDisp);

gpFit.hyp=hyp;
gpFit.optHyp = optHyp;
gpFit.xDisp = xDisp;
gpFit.yMu = yMu;
gpFit.yStd = sqrt(yVar);
gpFit.newPredFunc = @(X,Y,xDisp) gp(optHyp,infApx,meanFunc,covFunc,likFunc,X,Y,xDisp);