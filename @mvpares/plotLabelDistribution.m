function plotLabelDistribution(h,data,groups,xLabelStr,yLabelStr,titleStr,y_lim)
% Method for plotting distributions of predicted labels

distributionPlot(h,data,'groups',groups,'color',...
    [245/255 222/255 179/255],'showMM',0);
plotSpread(h,data,'distributionIdx',groups,...
    'distributionColors',[165/255 42/255 42/255]);
title(titleStr);
ylim(y_lim);
if ~isempty(xLabelStr), xlabel(xLabelStr); end
if ~isempty(yLabelStr), ylabel(yLabelStr); end

end