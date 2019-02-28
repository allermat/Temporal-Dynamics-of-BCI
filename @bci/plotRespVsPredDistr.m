function hFig = plotRespVsPredDistr(respDistr,predDistr)
% Method for plotting true vs predictied continuous distributions

nConds = size(respDistr,1);

sizeSubPlots = numSubplots(nConds);

hFig = figure();

for iCond = 1:nConds
    
    subplot(sizeSubPlots(1),sizeSubPlots(2),iCond);
    
    actRespF = respDistr.f_respA{iCond};
    actRespX = respDistr.x_respA{iCond};
    plot(actRespX,actRespF,...
        'LineWidth',1.5,'Color','r'); hold on;
    xMax = max(actRespX(actRespF > 1e-5));
    xMin = min(actRespX(actRespF > 1e-5));
    actRespF = respDistr.f_respV{iCond};
    actRespX = respDistr.x_respV{iCond};
    plot(actRespX,actRespF,...
        'LineWidth',1.5,'Color','b');
    xMax = max([xMax,max(actRespX(actRespF > 1e-5))]); 
    xMin = min([xMin,min(actRespX(actRespF > 1e-5))]);
    actPredF = predDistr.f_predA{iCond};
    actPredX = predDistr.x_predA{iCond};
    plot(actPredX,actPredF,...
        'LineWidth',1.5,'Color','r','LineStyle','--');
    xMax = max([xMax,max(actPredX(actPredF > 1e-5))]); 
    xMin = min([xMin,min(actPredX(actPredF > 1e-5))]);
    actPredF = predDistr.f_predV{iCond};
    actPredX = predDistr.x_predV{iCond};
    plot(actPredX,actPredF,...
        'LineWidth',1.5,'Color','b','LineStyle','--');
    xMax = max([xMax,max(actPredX(actPredF > 1e-5))]); 
    xMin = min([xMin,min(actPredX(actPredF > 1e-5))]);
    
%     xMin = min(actPredX(actPredF > 1e-5));
%     xMax = max(actPredX(actPredF > 1e-5));
    
    xlim([xMin,xMax]);
    
    if iCond == 1
        legend({'A resp','V resp','A pred','V pred'},'Location','NorthWestOutside');
    end
end

set(gcf,'Units','normalized','Position',[0,0,1,1],...
    'PaperPositionMode','auto');
suplabel(sprintf('Response vs predicted distributions'),'t',[0.08,0.08,0.84,0.86]);

end

