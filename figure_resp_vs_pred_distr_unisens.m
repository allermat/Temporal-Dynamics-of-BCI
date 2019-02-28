function hFig = figure_resp_vs_pred_distr_unisens(subID)

% Parsing input, checking matlab
p = inputParser;
addRequired(p,'subID',@ischar);
parse(p,subID);
subID = p.Results.subID;

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
    warning('analysis_behav_descriptives:missingData',...
        'Missing behavioural data, returning');
    hFig = [];
    return;
end

if strcmp(subID,'group')
    filePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_group'),...
        ['bci_simul_BEHAV_',subID,'.mat']);
else
    filePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),...
        ['bci_simul_BEHAV_',subID,'.mat']);
end
if exist(filePath,'file')
    bciSimulations = load(filePath);
    bciSimulations = bciSimulations.bciSimulations;
    predDistribution = cat(1,bciSimulations.conditions);
    predDistribution.relV = [bciSimulations.relV]';
    predDistribution.freq_predV = cat(1,bciSimulations.freq_predV);
    predDistribution.freq_predA = cat(1,bciSimulations.freq_predA);
else
    warning('analysis_behav_descriptives:missingData',...
        'Missing simulation data, returning');
    hFig = [];
    return;
end

% Taking respDistribution into unisensory and bisensory pieces
respDistribution.uniSens = false(size(respDistribution,1),1);
respDistribution.uniSens(isnan(respDistribution.locV) | ...
                         isnan(respDistribution.locA)) = true;
locLevels = unique(respDistribution.locV(~isnan(respDistribution.locV)));
visLocLevels = locLevels;
audLocLevels = locLevels';
visLocLevels = cat(1,NaN(1,size(visLocLevels,1)+1),repmat(visLocLevels,1,5))';
audLocLevels = cat(2,NaN(size(audLocLevels,2)+1,1),repmat(audLocLevels,5,1))';
visRelLevels = unique(respDistribution.relV(~isnan(respDistribution.relV)));
taskLevels = unique(respDistribution.task);
[nRows,nCols] = deal(numel(locLevels)+1); % We plot the
                                             % unisensories as well

for iRel = 1:numel(visRelLevels)
    hFig(iRel) = figure();
    iSubplot = 1;
    for iRow = 1:nRows
        for iCol = 1:nCols
            
            if iSubplot == 1
                % Skipping top left corner
                iSubplot = iSubplot + 1;
                continue;
            end
            
            subplot(nRows,nCols,iSubplot);
            if iRow == 1 && iCol ~= 1
                actRespDistr = respDistribution.distribution(...
                    respDistribution.locA == audLocLevels(iSubplot) &...
                    respDistribution.uniSens & ...
                    respDistribution.task == taskLevels(1),:);
                plot(locLevels,actRespDistr,...
                     'LineWidth',1.5,'Color','r');
            elseif iRow ~= 1 && iCol == 1
                actRespDistr = respDistribution.distribution(...
                    respDistribution.locV == visLocLevels(iSubplot) &...
                    respDistribution.relV == visRelLevels(iRel) &...
                    respDistribution.uniSens & ...
                    respDistribution.task == taskLevels(2),:);
                plot(locLevels,actRespDistr,...
                     'LineWidth',1.5,'Color','b');
            else
                actRespDistr = respDistribution.distribution(...
                    respDistribution.locV == visLocLevels(iSubplot) &...
                    respDistribution.locA == audLocLevels(iSubplot) &...
                    respDistribution.relV == visRelLevels(iRel) &...
                    respDistribution.task == taskLevels(1),:);
                plot(locLevels,actRespDistr,...
                     'LineWidth',1.5,'Color','r'); hold on;
                actRespDistr = respDistribution.distribution(...
                    respDistribution.locV == visLocLevels(iSubplot) &...
                    respDistribution.locA == audLocLevels(iSubplot) &...
                    respDistribution.relV == visRelLevels(iRel) &...
                    respDistribution.task == taskLevels(2),:);
                plot(locLevels,actRespDistr,...
                     'LineWidth',1.5,'Color','b');
                actPredDistr = predDistribution.freq_predA(...
                    predDistribution.locV == visLocLevels(iSubplot) &...
                    predDistribution.locA == audLocLevels(iSubplot) &...
                    predDistribution.relV == visRelLevels(iRel),:);
                plot(locLevels,actPredDistr,...
                     'LineWidth',1.5,'Color','r','LineStyle','--');
                actPredDistr = predDistribution.freq_predV(...
                    predDistribution.locV == visLocLevels(iSubplot) &...
                    predDistribution.locA == audLocLevels(iSubplot) &...
                    predDistribution.relV == visRelLevels(iRel),:);
                plot(locLevels,actPredDistr,...
                     'LineWidth',1.5,'Color','b','LineStyle','--');
            end
            
            ylim([0 1]);
            if iRow == 1
                set(gca,'XTick',[]);
                if iCol == 2
                    ylabel('No Vision');
                end
            elseif iRow == nRows && iCol ~= 1
                xlabel(sprintf('\n%.2f',audLocLevels(iSubplot)));
                set(gca,'XTick',locLevels);
            else
                set(gca,'XTick',[]);
            end
            
            if iCol == 1
                ylabel(sprintf('%.2f\n',visLocLevels(iSubplot)));
                if iRow == 2
                    title('No Audio','FontWeight','normal');
                end
            else
                set(gca,'YTickLabel',[]);
            end
            
            if iSubplot == numel(visLocLevels)
                legend({'A resp','V resp','A pred','V pred'},...
                       'Location','Northwest');
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
    suplabel(sprintf('Response vs predicted distributions, %s VR \n %s',relStr,subjStr),...
        't',[0.08,0.08,0.84,0.86]);
    suplabel('Visual location','y');
    suplabel('Auditory location','x');
end

end

