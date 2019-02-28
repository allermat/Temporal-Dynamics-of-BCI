function hFig = showAVmodelWeights(obj,effect,varargin)
% Method for plotting estimated AV model weigths
% 
% USAGE: 
%   hFig = showAVmodelWeights(obj,effect)
%   hFig = showAVmodelWeights(obj,effect,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%       effect (string): indicating which effects to plot. Possible values:
%           'int3way','int2way','main' respectively for 3-way interaction,
%           2-way interactions, main effects.
%   'Name'-Value arguments:
%       addStats (logical): wheter to add markers indicating significant
%           effects. Default: true.
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'.
%       plotError (logical): wheter to add error bars (takes effect only if
%           time courses are plotted (e.g. genTime is 'tr'). Default: true
%       smooth (logical): whether to smooth the time course of AV estimates
%       cust (struct): various settings to customize the plot's
%           appearence
%           Possible fields: supTitle, degLims
% OUTPUT:
%   hFig (scalar/vector): figure handles for plotted figures  

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validEffects = {'int3way','int2way','main','main_diff'};
validGenTimes = {'tr','tr_x_tr'};
addRequired(p,'obj');
addRequired(p,'effect',@(x)any(validatestring(x,validEffects)));
addParameter(p,'addStats',true,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'genTime','tr',@(x) any(validatestring(x,validGenTimes)));
addParameter(p,'plotError',true,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'smooth',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'cust',[],@isstruct);

parse(p,obj,effect,varargin{:});

obj = p.Results.obj;
effect = p.Results.effect;
addStats = p.Results.addStats;
genTime = p.Results.genTime;
plotError = p.Results.plotError;
smooth = p.Results.smooth;
cust = p.Results.cust;

weights = obj.getAVmodelWeights('genTime',genTime,'smooth',smooth,'unit','deg');
if isempty(weights)
    warning('mvpares:showAVmodelWeights:reqestedDataNotPresent',...
        'Couldn''t find AV model weigths in the dataset, returning.');
    return;
end
weightFnames = fieldnames(weights);
if addStats
    stats = obj.getStats('avWeights','genTime',genTime,'smooth',smooth);
    if isempty(stats)
        warning('mvpares:showAVmodelWeights:reqestedDataNotPresent',...
            ['Couldn''t find statistics for AV model weigths in the dataset, ',...
            'plotting without statistics.']);
        addStats = false;
    else
        statFields = fieldnames(stats);
        stats = struct2cell(stats);
    end
else
    stats = [];
end

hFig = [];

time = obj.getTrTimePoints;
sTime = obj.plotTimeWin(1);
eTime = obj.plotTimeWin(2);
if isfield(cust,'degLims') && ~isempty(cust.degLims)
    degLims = cust.degLims;
else
    degLims = [-20,150];
end

switch effect
    case 'int3way'
        hActFig = figure();
        hFig = [hFig,hActFig];
        det.hFig = hActFig;
        det.dataUnit = 'W_A_V (degree)';
        if addStats
            if all(cellfun(@isempty,regexp(statFields,'h_VR_X_Task_X_Disp','once')))
                warning('mvpares:showAVmodelWeights:reqestedDataNotPresent',...
                        ['Couldn''t find statistics for AV model weigths for ',...
                         'this particular effect in the dataset, plotting ',...
                         'without statistics.']);
                addStats = false;
            end
        end
        switch genTime
            case 'tr'
                data = cat(2,...
                    weights.r_d_a_wav,weights.r_d_v_wav,...
                    weights.r_D_a_wav,weights.r_D_v_wav,...
                    weights.R_d_a_wav,weights.R_d_v_wav,...
                    weights.R_D_a_wav,weights.R_D_v_wav);
                error = cat(2,...
                    weights.r_d_a_ci,weights.r_d_v_ci,...
                    weights.r_D_a_ci,weights.r_D_v_ci,...
                    weights.R_d_a_ci,weights.R_d_v_ci,...
                    weights.R_D_a_ci,weights.R_D_v_ci);
                if addStats
                    idx = ~cellfun(@isempty,regexp(statFields,...
                                                   'h_VR_X_Task_X_Disp','once'));
                    if any(idx)
                        actH = stats{idx};
                        clusterIdx = mvpares.findClusters(actH);
                        det.addRect = arrayfun(@(x) time(x),clusterIdx);
                    else
                        warning('mvpares:showAVmodelWeights:reqestedDataNotPresent',...
                                ['Couldn''t find statistics for this effect ',...
                                 'in the dataset, plotting without statistics.']);
                        det.addRect = [];
                    end
                end
                det.lineProp = cat(1,...
                    {'Color','k','LineWidth',1.5,'LineStyle','-'},...
                    {'Color','k','LineWidth',1.5,'LineStyle','--'},...
                    {'Color','k','LineWidth',1.5,'LineStyle',':'},...
                    {'Color','k','LineWidth',1.5,'LineStyle','-.'},...
                    {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-'},...
                    {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','--'},...
                    {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle',':'},...
                    {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-.'});
                det.legend = {'VR-, D-, A','VR-, D-, V','VR-, D+, A','VR-, D+, V','VR+, D-, A',...
                    'VR+, D-, V','VR+, D+, A','VR+, D+, V','Location','SouthEastOutside'};
                det.plotError = plotError;
                det.title = 'W_A_V, 3-way interaction';
                det.yLim = degLims;
                det.xLim = [sTime,eTime];
                det.addLine = {[[sTime eTime]' [sTime eTime]'],[[90 90]' [0 0]'],...
                        'Color',[0.5 0.5 0.5]};
                hFig = mvpares.plotTimeSeriesData(time,data,error,det);
                set(hActFig,'Units','normalized','Position',[1/4,1/4,1/2,1/2],...
                    'PaperPositionMode','auto');
            case 'tr_x_tr'
                condsOfInterest = {'R_D_a','R_D_v','R_d_a','R_d_v',...
                                   'r_D_a','r_D_v','r_d_a','r_d_v'};
                titles = {'VR+, D+, A','VR+, D+, V','VR+, D-, A','VR+, D-, V',...
                          'VR-, D+, A','VR-, D+, V','VR-, D-, A','VR-, D-, V'};
                nColFig = 4;
                nRowFig = ceil(numel(condsOfInterest)/nColFig);
                for i = 1:numel(condsOfInterest)
                    data = weights.(strcat(condsOfInterest{i},'_wav'));
                    det.hAxes = subplot(nRowFig,nColFig,i);
                    det.title = titles{i};
                    det.cLim = degLims;
                    det.xLim = [sTime,eTime];
                    mvpares.plotTimeByTimeMatrix(time,data,det);
                end
                set(hFig,'Units','normalized','Position',[1/8,1/4,3/4,0.6],...
                    'PaperPositionMode','auto');
                temp = findobj(hFig,'Type','axes');
                hAxes = temp(1);
                posOutAxes = get(hAxes,'OuterPosition');
                posOutCbar = [posOutAxes(1)+(posOutAxes(3)),posOutAxes(2),...
                    (posOutAxes(3)*0.3),posOutAxes(4)];
                hCbar = colorbar('OuterPosition',posOutCbar);
                ylabel(hCbar,'W_A_V (degree)');
                if isfield(cust,'supTitle')
                    suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
                else
                    suplabel('W_A_V 3-way interaction','t',[0.08,0.08,0.84,0.88]);
                end
        end
        
    case 'int2way'
        hActFig = figure();
        hFig = [hFig,hActFig];
        det.hFig = hActFig;
        det.dataUnit = 'W_A_V (degree)';
        condsOfInterest = {...
                    'r_d','r_D','R_d','R_D';...
                    'r_a','r_v','R_a','R_v';...
                    'd_a','d_v','D_a','D_v'};
        isCondPresent = true(size(condsOfInterest));
        for i = 1:numel(condsOfInterest)
            isCondPresent(i) = any(~cellfun(@isempty, ...
                                       regexp(weightFnames,['^',condsOfInterest{i},'_wav'])));
        end
        % The length of this vector will be size(isCondPresent,1)
        isEffectPresent = all(isCondPresent,2);
        % condsOfInterest = condsOfInterest(isCondPresent);
        % The array needs to be rashaped to an N x 4 matrix
        % condsOfInterest = reshape(condsOfInterest,[],4);
        if addStats
            if all(cellfun(@isempty,...
                           regexp(statFields,'h_(VR_X_Disp|VR_X_Task|Disp_X_Task)','once')))
                warning('mvpares:showAVmodelWeights:reqestedDataNotPresent',...
                        ['Couldn''t find statistics for AV model weigths for ',...
                         'this particular effect in the dataset, plotting ',...
                         'without statistics.']);
                addStats = false;
            end
        end
        switch genTime
            case 'tr'
                titles = {'VR X Disparity','VR X Task','Disparity X Task'};
                effectNames = {'VR_X_Disp','VR_X_Task','Disp_X_Task'};
                legends = {...
                    'VR-, D-','VR-, D+','VR+, D-','VR+, D+','Location','SouthEastOutside';...
                    'VR-, A','VR-, V','VR+, A','VR+, V','Location','SouthEastOutside';...
                    'D-, A','D-, V','D+, A','D+, V','Location','SouthEastOutside'};
                % titles = titles(isEffectPresent);
                % legends = legends(isEffectPresent,:);
                for i = 1:numel(titles)
                    if ~isEffectPresent(i), continue; end
                    actConds = condsOfInterest(i,:);
                    actWavs = strcat(actConds,'_wav');
                    actCis = strcat(actConds,'_ci');
                    data = cat(2,...
                        weights.(actWavs{1}),weights.(actWavs{2}),...
                        weights.(actWavs{3}),weights.(actWavs{4}));
                    error = cat(2,...
                        weights.(actCis{1}),weights.(actCis{2}),...
                        weights.(actCis{3}),weights.(actCis{4}));
                    if addStats
                        idx = ~cellfun(@isempty,regexp(statFields,...
                                                       ['h_',effectNames{i}],'once'));
                        if any(idx)
                            actH = stats{idx};
                            clusterIdx = mvpares.findClusters(actH);
                            det.addRect = arrayfun(@(x) time(x),clusterIdx);
                        else
                            det.addRect = [];
                        end
                    end
                    det.lineProp = cat(1,...
                        {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-'},...
                        {'Color','k','LineWidth',1.5,'LineStyle','-'},...
                        {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','--'},...
                        {'Color','k','LineWidth',1.5,'LineStyle','--'});
                    det.hAxes = subplot(numel(titles),1,i);
                    det.legend = legends(i,:);
                    det.plotError = plotError;
                    det.title = titles{i};
                    det.yLim = degLims;
                    det.xLim = [sTime,eTime];
                    det.addLine = {[[sTime eTime]' [sTime eTime]'],[[90 90]' [0 0]'],...
                        'Color',[0.5 0.5 0.5]};
                    mvpares.plotTimeSeriesData(time,data,error,det);
                end
                % Setting figure size
                set(hActFig,'Units','normalized','Position',[1/3,0,1/3,1],...
                    'PaperPositionMode','auto');
                if isfield(cust,'supTitle')
                    suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
                else
                    suplabel('W_A_V 2-way interactions','t',[0.08,0.08,0.84,0.88]);
                end
            case 'tr_x_tr'
                titles = {'VR-, D-','VR-, D+','VR+, D-','VR+, D+',...
                          'VR-, A','VR-, V','VR+, A','VR+, V',...
                          'D-, A','D-, V','D+, A','D+, V'};
                nColFig = 4;
                nRowFig = ceil(numel(condsOfInterest)/nColFig);
                for i = 1:numel(condsOfInterest)
                    data = weights.(strcat(condsOfInterest{i},'_wav'));
                    det.hAxes = subplot(nRowFig,nColFig,i);
                    det.title = titles{i};
                    det.cLim = degLims;
                    det.xLim = [sTime,eTime];
                    mvpares.plotTimeByTimeMatrix(time,data,det);
                end
                set(hFig,'Units','normalized','Position',[1/8,0,3/4,0.9],...
                    'PaperPositionMode','auto');
                temp = findobj(hFig,'Type','axes');
                hAxes = temp(1);
                posOutAxes = get(hAxes,'OuterPosition');
                posOutCbar = [posOutAxes(1)+(posOutAxes(3)),posOutAxes(2),...
                    (posOutAxes(3)*0.3),posOutAxes(4)];
                hCbar = colorbar('OuterPosition',posOutCbar);
                ylabel(hCbar,'W_A_V (degree)');
                if isfield(cust,'supTitle')
                    suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
                else
                    suplabel('W_A_V 2-way interactions','t',[0.08,0.08,0.84,0.88]);
                end
        end
        
    case 'main'
        hActFig = figure();
        hFig = [hFig,hActFig];
        det.hFig = hActFig;
        det.dataUnit = 'W_A_V (degree)';
        condsOfInterest = {'r','R';'a','v';'d','D'};
        if addStats
            if all(cellfun(@isempty,regexp(statFields,'h_(Disp|Task|VR)','once')))
                warning('mvpares:showAVmodelWeights:reqestedDataNotPresent',...
                    ['Couldn''t find statistics for AV model weigths for ',...
                    'this particular effect in the dataset, plotting ',...
                    'without statistics.']);
                addStats = false;
            end
        end
        % Checking which conditions are present in the dataset
        isCondPresent = true(size(condsOfInterest));
        for i = 1:numel(condsOfInterest)
            isCondPresent(i) = any(~cellfun(@isempty, ...
                                       regexp(weightFnames,['^',condsOfInterest{i},'_wav'])));
        end
        % The length of this vector will be size(isCondPresent,1)
        isEffectPresent = all(isCondPresent,2);
        % condsOfInterest = condsOfInterest(isCondPresent);
        % The array needs to be rashaped to an N x 2 matrix
        % condsOfInterest = reshape(condsOfInterest,[],2);
        switch genTime
          case 'tr'
            % Defining various variables related to the effects
            titles = {'Visual reliability','Task','Disparity'};
            effectNames = {'VR','Task','Disp'};
            legends = {...
                'VR-','VR+','Location','SouthEastOutside';...
                'A','V','Location','SouthEastOutside';...
                'D-','D+','Location','SouthEastOutside'};
            lineProps = {...
                cat(1,...
                    {'Color','k','LineWidth',1.5,'LineStyle','-'},...
                    {'Color','k','LineWidth',1.5,'LineStyle','--'}),...
                cat(1,...
                    {'Color',[0.5 0.5 0.5],'LineWidth',1.5},...
                    {'Color','k','LineWidth',1.5}),...
                cat(1,...
                    {'Color','k','LineWidth',1.5,'LineStyle','-'},...
                    {'Color','k','LineWidth',1.5,'LineStyle','--'})};
            % titles = titles(isEffectPresent);
            % effectNames = effectNames(isEffectPresent);
            % legends = legends(isEffectPresent,:);
            % lineProps = lineProps(isEffectPresent);
            for i = 1:numel(effectNames)                
                % Simply skip those effects which are not present
                if ~isEffectPresent(i), continue; end
                actConds = condsOfInterest(i,:);
                actWavs = strcat(actConds,'_wav');
                actCis = strcat(actConds,'_ci');
                data = cat(2,weights.(actWavs{1}),weights.(actWavs{2}));
                error = cat(2,weights.(actCis{1}),weights.(actCis{2}));
                if addStats
                    idx = ~cellfun(@isempty,regexp(statFields,...
                                                   ['h_',effectNames{i}],'once'));
                    if any(idx)
                        actH = stats{idx};
                        clusterIdx = mvpares.findClusters(actH);
                        det.addRect = arrayfun(@(x) time(x),clusterIdx);
                    else
                        det.addRect = [];
                    end
                end
                det.lineProp = lineProps{i};
                det.hAxes = subplot(numel(effectNames),1,i);
                det.legend = legends(i,:);
                det.plotError = plotError;
                det.title = titles{i};
                det.yLim = degLims;
                det.xLim = [sTime,eTime];
                det.addLine = {[[sTime eTime]' [sTime eTime]'],[[90 90]' [0 0]'],...
                               'Color',[0.5 0.5 0.5]};
                mvpares.plotTimeSeriesData(time,data,error,det);
            end
            % Setting figure size
            set(hActFig,'Units','normalized','Position',[1/3,0,1/3,1],...
                        'PaperPositionMode','auto');
            if isfield(cust,'supTitle')
                suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
            else
                suplabel('W_A_V main effects','t',[0.08,0.08,0.84,0.88]);
            end
          case 'tr_x_tr'
            titles = {'VR-','VR+','A','V','D-','D+'};
            titles = titles(isCondPresent);
            [nRowFig,nColFig,iPlot] = mvpares.segmentFigArea(floor(numel(titles))/2,2,1/4);
            for i = 1:numel(condsOfInterest)
                data = weights.(strcat(condsOfInterest{i},'_wav'));
                det.hAxes = subplot(nRowFig,nColFig,iPlot{i});
                det.title = titles{i};
                det.cLim = degLims;
                det.xLim = [sTime,eTime];
                mvpares.plotTimeByTimeMatrix(time,data,det);
            end
            set(hFig,'Units','normalized','Position',[1/4,0,0.37,0.9],...
                     'PaperPositionMode','auto');
            temp = findobj(hFig,'Type','axes');
            hAxes = temp(1);
            posOutAxes = get(hAxes,'OuterPosition');
            posOutCbar = [posOutAxes(1)+(posOutAxes(3)),posOutAxes(2),...
                          (posOutAxes(3)*0.3),posOutAxes(4)];
            hCbar = colorbar('OuterPosition',posOutCbar);
            ylabel(hCbar,'W_A_V (degree)');
            if isfield(cust,'supTitle')
                suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
            else
                suplabel('W_A_V main effects','t',[0.08,0.08,0.84,0.88]);
            end
        end
            
    case 'main_diff'
        hActFig = figure();
        hFig = [hFig,hActFig];
        det.hFig = hActFig;
        det.dataUnit = 'W_A_V (degree)';
        condsOfInterest = {'R','r';'v','a';'D','d'};
        switch genTime
            case 'tr'
%                 titles = {'Visual Reliability','Task','Disparity'};
%                 legends = {...
%                     'VR-','VR+','Location','SouthEastOutside';...
%                     'A','V','Location','SouthEastOutside';...
%                     'D-','D+','Location','SouthEastOutside'};
%                 lineProps = {...
%                     cat(1,...
%                     {'Color','k','LineWidth',1.5,'LineStyle','-'},...
%                     {'Color','k','LineWidth',1.5,'LineStyle','--'}),...
%                     cat(1,...
%                     {'Color',[0.5 0.5 0.5],'LineWidth',1.5},...
%                     {'Color','k','LineWidth',1.5}),...
%                     cat(1,...
%                     {'Color','k','LineWidth',1.5,'LineStyle','-'},...
%                     {'Color','k','LineWidth',1.5,'LineStyle','--'})};
%                 for i = 1:numel(titles)
%                     actConds = condsOfInterest(i,:);
%                     actWavs = strcat(actConds,'_wav');
%                     actCis = strcat(actConds,'_ci');
%                     data = cat(2,weights.(actWavs{1}),weights.(actWavs{2}));
%                     error = cat(2,weights.(actCis{1}),weights.(actCis{2}));
%                     det.lineProp = lineProps{i};
%                     det.hAxes = subplot(3,1,i);
%                     det.legend = legends(i,:);
%                     det.plotError = plotError;
%                     det.title = titles{i};
%                     det.yLim = degLims;
%                     det.xLim = [sTime,eTime];
%                     det.addLine = {[[sTime eTime]' [sTime eTime]'],[[90 90]' [0 0]'],...
%                         'Color',[0.5 0.5 0.5]};
%                     mvpares.plotTimeSeriesData(time,data,error,det);
%                 end
%                 % Setting figure size
%                 set(hActFig,'Units','normalized','Position',[1/3,0,1/3,1],...
%                     'PaperPositionMode','auto');
%                 suplabel('W_A_V main effects','t',[0.08,0.08,0.84,0.88]);
            case 'tr_x_tr'
                titles = {'(VR+) - (VR-)','V - A','(D+) - (D-)'};
                [nRowFig,nColFig,iPlot] = mvpares.segmentFigArea(3,1,1/4);
                for i = 1:numel(titles)
                    data1 = weights.(strcat(condsOfInterest{i,1},'_wav'));
                    data2 = weights.(strcat(condsOfInterest{i,2},'_wav'));
                    data = data1-data2;
                    det.hAxes = subplot(nRowFig,nColFig,iPlot{i});
                    det.title = titles{i};
                    det.cLim = [-80,80];
                    det.xLim = [sTime,eTime];
                    mvpares.plotTimeByTimeMatrix(time,data,det);
                end
                set(hActFig,'Units','normalized','Position',[1/3,0,0.2,0.9],...
                    'PaperPositionMode','auto');
                temp = findobj(hFig,'Type','axes');
                hAxes = temp(1);
                posOutAxes = get(hAxes,'OuterPosition');
                posOutCbar = [posOutAxes(1)+(posOutAxes(3)),posOutAxes(2),...
                    (posOutAxes(3)*0.3),posOutAxes(4)];
                hCbar = colorbar('OuterPosition',posOutCbar);
                ylabel(hCbar,'Difference (degree)');
                if isfield(cust,'supTitle')
                    suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
                else
                    suplabel('W_A_V differences, main effects','t',[0.08,0.08,0.84,0.88]);
                end
        end
end

end