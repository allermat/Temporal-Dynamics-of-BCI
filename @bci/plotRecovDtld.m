function hFig = plotRecovDtld(respDistr,predDistr1,predDistr2,varargin)
% Method for plotting true vs predictied continuous distributions

p = inputParser;
addRequired(p,'respDistr',@(x) validateattributes(x,{'table'},{'nonempty'}));
addRequired(p,'predDistr1',@(x) validateattributes(x,{'table'},{'nonempty'}));
addRequired(p,'predDistr2',@(x) validateattributes(x,{'table'},{'nonempty'}));
addOptional(p,'nPlots',size(respDistr,1),@(x) validateattributes(x,{'numeric'},...
    {'scalar','integer','positive','<=',size(respDistr,1)}));
addOptional(p,'titleStr','',@(x) validateattributes(x,{'char'},{'nonempty'}));

parse(p,respDistr,predDistr1,predDistr2,varargin{:});

respDistr = p.Results.respDistr;
predDistr1 = p.Results.predDistr1;
predDistr2 = p.Results.predDistr2;
nPlots = p.Results.nPlots;
titleStr = p.Results.titleStr;

for iPlot = 1:nPlots
    
    hFig(iPlot) = figure();
        
    actRespF = respDistr.f_respA{iPlot};
    actRespX = respDistr.x_respA{iPlot};
    plot(actRespX,actRespF,...
        'LineWidth',1.5,'Color','r'); hold on;
    xMax = max(actRespX(actRespF > 1e-5));
    xMin = min(actRespX(actRespF > 1e-5));
    actRespF = respDistr.f_respV{iPlot};
    actRespX = respDistr.x_respV{iPlot};
    plot(actRespX,actRespF,...
        'LineWidth',1.5,'Color','b');
    xMax = max([xMax,max(actRespX(actRespF > 1e-5))]); 
    xMin = min([xMin,min(actRespX(actRespF > 1e-5))]);
    actPredF = predDistr1.f_predA{iPlot};
    actPredX = predDistr1.x_predA{iPlot};
    plot(actPredX,actPredF,...
        'LineWidth',1.5,'Color','r','LineStyle','--');
    xMax = max([xMax,max(actPredX(actPredF > 1e-5))]); 
    xMin = min([xMin,min(actPredX(actPredF > 1e-5))]);
    actPredF = predDistr1.f_predV{iPlot};
    actPredX = predDistr1.x_predV{iPlot};
    plot(actPredX,actPredF,...
        'LineWidth',1.5,'Color','b','LineStyle','--');
    xMax = max([xMax,max(actPredX(actPredF > 1e-5))]); 
    xMin = min([xMin,min(actPredX(actPredF > 1e-5))]);
    actPredF = predDistr2.f_predA{iPlot};
    actPredX = predDistr2.x_predA{iPlot};
    plot(actPredX,actPredF,...
        'LineWidth',1.5,'Color','r','LineStyle','-.');
    xMax = max([xMax,max(actPredX(actPredF > 1e-5))]); 
    xMin = min([xMin,min(actPredX(actPredF > 1e-5))]);
    actPredF = predDistr2.f_predV{iPlot};
    actPredX = predDistr2.x_predV{iPlot};
    plot(actPredX,actPredF,...
        'LineWidth',1.5,'Color','b','LineStyle','-.');
    xMax = max([xMax,max(actPredX(actPredF > 1e-5))]); 
    xMin = min([xMin,min(actPredX(actPredF > 1e-5))]);
%     xMin = min(actPredX(actPredF > 1e-5));
%     xMax = max(actPredX(actPredF > 1e-5));
    
    xlim([xMin,xMax]);
    
    if iPlot == 1
        legend({'A Resp_p_s_e_u','V Resp_p_s_e_u','sA\_resp|Par_T_R_U_E','sV\_resp|Par_T_R_U_E','sA\_resp|Par_E_S_T','sV\_resp|Par_E_S_T'});
    end
    title(sprintf('sV: %.2f, sA: %.2f \n%s',respDistr.locV(iPlot),respDistr.locA(iPlot),titleStr));
end

end

