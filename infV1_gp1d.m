function gpFit = infV1_gp1d(X,Y)

gpFit = struct();

hyp = struct();
likFunc = {@likGauss};
hyp.lik = log(std(Y));
covFunc = {@covSEiso};
hyp.cov = [log(std(X));0];
meanFunc = {@meanConst};
hyp.mean = 0;
infFunc = @infGaussLik;

xu = linspace(prctile(X,0.1),prctile(X,99.9),30)'; 
covFunc = {'apxSparse',covFunc,xu};      % inducing points
infApx  = @(varargin) infFunc(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
hyp.xu = xu;

optHyp = minimize(hyp,@gp,-150,infApx,meanFunc,covFunc,likFunc,X,Y);

xDisp = linspace(prctile(X,0.01),prctile(X,99.99),1e3)';
[~,~,yMu,yVar] = gp(optHyp,infApx,meanFunc,covFunc,likFunc,X,Y,xDisp);

gpFit.hyp=hyp;
gpFit.optHyp = optHyp;
gpFit.xDisp = xDisp;
gpFit.yMu = yMu;
gpFit.yStd = sqrt(yVar);