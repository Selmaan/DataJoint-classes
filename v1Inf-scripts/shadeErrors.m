function hL = shadeErrors(y, col, x, ste_flag)

if nargin<3
    x = 1:size(y,2);
end

if nargin<4
    ste_flag = true;
end

mVec = mean(y);
if ste_flag
    sVec = std(y)./sqrt(size(y,1));
else
    sVec = std(y);
end

mVec = mVec(:); sVec = sVec(:); x = x(:);
patch([x;flipud(x)],[mVec+sVec; flipud(mVec-sVec)],col,...
    'FaceAlpha',0.25,'EdgeColor','none'),
hL = plot(x,mVec,'color',col,'linewidth',2);