function gpFit = infV1_gp2d(X,Y)

gpFit = struct();

hyp = struct();
likFunc = {@likGauss};
hyp.lik = log(std(Y));
% covFunc = {@covADD,{[1,2],'covSEiso'}};
% hyp.cov = [log(std(X(:,1)));0;log(std(X(:,2)));0;0;0];
covFunc_dist = {@covMask {[1 0], @covSEiso}};
covFunc_corr = {@covMask {[0 1], @covSEiso}};
covFunc = {@covSum {covFunc_dist, covFunc_corr}};
hyp.cov = [log(std(X(:,1)));0;log(std(X(:,2)));0];

meanFunc = {@meanConst};
hyp.mean = 0;
infFunc = @infGaussLik;

xu_1 = linspace(prctile(X(:,1),0.1),prctile(X(:,1),99.9),10)'; 
xu_2 = linspace(prctile(X(:,2),0.1),prctile(X(:,2),99.9),10)'; 
[xu_1,xu_2] = ndgrid(xu_1,xu_2);
xu = [xu_1(:), xu_2(:)];

covFunc = {'apxSparse',covFunc,xu};      % inducing points
infApx  = @(varargin) infFunc(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
hyp.xu = xu;

optHyp = minimize(hyp,@gp,-150,infApx,meanFunc,covFunc,likFunc,X,Y);

[disp_1, disp_2] = ndgrid(linspace(prctile(X(:,1),0.1),prctile(X(:,1),99.9),1e2)',...
    linspace(prctile(X(:,2),0.1),prctile(X(:,2),99.9),1e2)');
xDisp = [disp_1(:), disp_2(:)];
[~,~,yMu,yVar] = gp(optHyp,infApx,meanFunc,covFunc,likFunc,X,Y,xDisp);

gpFit.hyp=hyp;
gpFit.optHyp = optHyp;
gpFit.xDisp = xDisp;
gpFit.yMu = yMu;
gpFit.yStd = sqrt(yVar);