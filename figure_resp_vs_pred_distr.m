function hFig = figure_resp_vs_pred_distr(subID,varargin)

% Parsing input
p = inputParser;

validModels = {'bci','fus','segA','segV','taskRel','null'};

addRequired(p,'subID',@ischar);
addOptional(p,'model','bci',@(x) any(validatestring(x,validModels)));

parse(p,subID,varargin{:});

subID = p.Results.subID;
model = p.Results.model;

% Checking if there is existing processed behavioural data
expStage = 'final';

if strcmp(subID,'group')
    filePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_group'),...
        ['descriptives_BEHAV_',subID,'.mat']);
else
    filePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),...
        ['descriptives_BEHAV_',subID,'.mat']);
end
if exist(filePath,'file')
    respDistribution = load(filePath);
    respDistribution = respDistribution.respDistribution;
else
    warning('figure_resp_vs_pred_distr:missingData',...
        'Missing behavioural data, returning');
    hFig = [];
    return;
end


if strcmp(subID,'group')
    fileMatchStr = ['bci_simul_BEHAV_',subID,'.mat'];
    saveDf = cd(DEC_2_setupdir(expStage,'anal_behav_group',subID));
else
    fileMatchStr = ['bci_simul_BEHAV_',subID,'_[0-9]+.mat'];
    saveDf = cd(DEC_2_setupdir(expStage,'anal_behav_sub',subID));
end
fileList = cellstr(ls);
matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));

if sum(matchID) == 0
    warning('figure_resp_vs_pred_distr:missingData',...
            'Missing simulation data, returning');
    cd(saveDf);
    hFig = [];
    return;
elseif sum(matchID) > 1
    cd(saveDf);
    warning('figure_resp_vs_pred_distr:missingData',...
            'Multiple simulation data files, returning');
    cd(saveDf);
    hFig = [];
    return;
else
    fileName = fileList(matchID);
    bciSimulations = load(fileName{:});
    if strcmp(subID,'group')
        bciSimulations = bciSimulations.bciSimulations;
    else
        models = bciSimulations.models;
        if any(ismember(models,{model}))
            bciSimulations = bciSimulations.mdlEval{ismember(models,{model})};
        else
            warning('figure_resp_vs_pred_distr:invalidInput',...
                    ['The specified model is not part of the dataset, ' ...
                     'returning. ']); 
            hFig = [];
            return;
        end
    end
end
cd(saveDf);

predDistribution = cat(1,bciSimulations.conditions);
predDistribution.relV = [bciSimulations.relV]';
predDistribution.freq_predV = cat(1,bciSimulations.freq_predV);
predDistribution.freq_predA = cat(1,bciSimulations.freq_predA);

visLocLevels = unique(respDistribution.locV(~isnan(respDistribution.locV)));
audLocLevels = visLocLevels;
visRelLevels = unique(respDistribution.relV(~isnan(respDistribution.relV)));
taskLevels = unique(respDistribution.task);
[nRows,nCols] = deal(numel(visLocLevels));
tol = eps('single');
for iRel = 1:numel(visRelLevels)
    hFig(iRel) = figure();
    iSubplot = 1;
    for iRow = 1:nRows
        for iCol = 1:nCols
            
            subplot(nRows,nCols,iSubplot);
            actRespDistr = respDistribution.distribution(...
                ismembertol(respDistribution.locV,visLocLevels(iRow),tol) &...
                ismembertol(respDistribution.locA,audLocLevels(iCol),tol) &...
                respDistribution.relV == visRelLevels(iRel) &...
                respDistribution.task == taskLevels(1),:);
            plot(visLocLevels,actRespDistr,...
                'LineWidth',1.5,'Color','r'); hold on;
            actRespDistr = respDistribution.distribution(...
                ismembertol(respDistribution.locV,visLocLevels(iRow),tol) &...
                ismembertol(respDistribution.locA,audLocLevels(iCol),tol) &...
                respDistribution.relV == visRelLevels(iRel) &...
                respDistribution.task == taskLevels(2),:);
            plot(visLocLevels,actRespDistr,...
                'LineWidth',1.5,'Color','b');
            actPredDistr = predDistribution.freq_predA(...
                ismembertol(predDistribution.locV,visLocLevels(iRow),tol) &...
                ismembertol(predDistribution.locA,audLocLevels(iCol),tol) &...
                predDistribution.relV == visRelLevels(iRel),:);
            plot(visLocLevels,actPredDistr,...
                'LineWidth',1.5,'Color','r','LineStyle','--');
            actPredDistr = predDistribution.freq_predV(...
                ismembertol(predDistribution.locV,visLocLevels(iRow),tol) &...
                ismembertol(predDistribution.locA,audLocLevels(iCol),tol) &...
                predDistribution.relV == visRelLevels(iRel),:);
            plot(visLocLevels,actPredDistr,...
                'LineWidth',1.5,'Color','b','LineStyle','--');
            
            ylim([0 1]);
            if iRow == 1
                set(gca,'XTick',[]);
            elseif iRow == nRows
                xlabel(sprintf('\n%.2f',audLocLevels(iCol)));
                set(gca,'XTick',audLocLevels);
            else
                set(gca,'XTick',[]);
            end
            
            if iCol == 1
                ylabel(sprintf('%.2f\n',visLocLevels(iRow)));
            else
                set(gca,'YTickLabel',[]);
            end
            
            if iSubplot == 1
                legend({'A resp','V resp','A pred','V pred'},'Location','NorthEast');
            end
            
            iSubplot = iSubplot + 1;
            
        end
    end
    
    switch iRel
        case 1
            relStr = 'high';
        case 2
            relStr = 'low';
    end
    
    if strcmp(subID,'group')
        subjStr = subID;
    else
        subjStr = ['subj: ',subID];
    end
    
    set(gcf,'Units','normalized','Position',[0.125,0.125,0.75,0.75],...
        'PaperPositionMode','auto');
    suplabel(sprintf('Response vs predicted distributions (%s), %s VR \n %s',model,relStr,subjStr),...
        't',[0.08,0.08,0.84,0.86]);
    suplabel('Visual location','y');
    suplabel('Auditory location','x');
end

end

