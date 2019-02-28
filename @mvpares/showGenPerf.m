function hFig = showGenPerf(obj,varargin)
% Method for plotting generalization prediction performance estimates
%
% USAGE:
%   hFig = showGenPerf(obj)
%   hFig = showGenPerf(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'.
%       perfIdx (string): performance index to be plotted.
%           Possible values: 
%           b - regression slope (beta),
%           r - correlation coeff,
%           rf - Fisher-transformed correlation coeff,
%           r2 - squared correlation coeff, 
%           acc - accuracy.
%           Default: r2 for regression, acc for classification
%       cvEval (string): cross-validaiton evaluation method.
%           Possible values: 'avgFolds','poolFolds'. Default:
%           'avgFolds'.
%       plotError (logical): wheter to add error bars (takes effect only if
%           time courses are plotted (e.g. genTime is 'tr'). Default: true
%       in (structure):
%
% OUTPUT:
%   hFig (scalar): figure handle for plotted figure
% 
% Copyright(C) 2016, Mate Aller

% Parsing input, checking matlab
p = inputParser;

validGenTimes = {'tr','tr_x_tr'};
validPerfIndices = {'r2','r','rf','b','acc'};
validCvEval = {'avgFolds','poolFolds'};

addRequired(p,'obj');
addParameter(p,'genTime','tr',@(x) any(validatestring(x,validGenTimes)));
addParameter(p,'perfIdx','',@(x) any(validatestring(x,validPerfIndices)));
addParameter(p,'cvEval','poolFolds',@(x) any(validatestring(x, ...
                                                  validCvEval)));
addParameter(p,'plotError',true,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'addStats',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
addParameter(p,'in',[],@isstruct);

parse(p,obj,varargin{:});

obj = p.Results.obj;
genTime = p.Results.genTime;
perfIdx = p.Results.perfIdx;
cvEval = p.Results.cvEval;
plotError = p.Results.plotError;
addStats = p.Results.addStats;
in = p.Results.in;

% Loading necessary data
svmType = obj.info.svm_type;
trCond = obj.info.tr_cond;
genCond = obj.info.gen_cond;
% Setting default performance index
if ismember(svmType,{'cc','nc'})
    if ismember(perfIdx,{'r','rf','r2','b'})
        warning('Changing accuracy index to acc!');
        perfIdx = 'acc';
    elseif ismember(perfIdx,{''})
        perfIdx = 'acc';
    end
elseif ismember(svmType,{'er','nr'})
    if ismember(perfIdx,{'acc'})
        warning('Changing accuracy index to r2!');
        perfIdx = 'rf';
    elseif ismember(perfIdx,{''})
        perfIdx = 'rf';
    end
end
% Loading requested performance estimate
if isprop(obj.data,'gen_perfEstimates')
    perfStruct = obj.data.gen_perfEstimates;
    fieldName = [perfIdx,'_',cvEval];
    if isfield(perfStruct,fieldName)
        perf = perfStruct.(fieldName);
    else
        hFig = [];
        warning('mvpares:showGenPerf:reqestedDataNotPresent',...
            ['The requested prediction performance estimate is not ',...
            'present in the dataset.']);
        return;
    end
else
    hFig = [];
    warning('mvpares:showGenPerf:reqestedDataNotPresent',...
        ['There are no prediciton performance estimates available ',...
        'in the dataset.']);
    return;
end

if addStats
    stats = obj.getStats('genPerf','genTime',genTime,'smooth',false);
    if isempty(stats)
        warning('mvpares:showGenPerf:reqestedDataNotPresent',...
            ['Couldn''t find statistics generalization performance in the dataset, ',...
            'plotting without statistics.']);
        addStats = false;
    else
        statFields = fieldnames(stats);
        stats = struct2cell(stats);
    end
else
    stats = [];
end
% Selecting time points
tr_timePoints = obj.getTrTimePoints;
gen_timePoints = tr_timePoints;

time = obj.getTrTimePoints;
sTime = obj.plotTimeWin(1);
eTime = obj.plotTimeWin(2);


if isfield(in,'lineProp'), lineProp = in.lineProp; else lineProp = {}; end
switch perfIdx
    case 'b', perfIdxStr = 'beta';
    case 'r', perfIdxStr = 'R';
    case 'rf', perfIdxStr = 'R (Fisher transformed)';
    case 'r2', perfIdxStr = 'R^2';
    case 'acc', perfIdxStr = 'Acc (%)';
end

% Plotting figure
if strcmp(genTime,'tr')
    
    if ~isvector(perf)
        perf = diag(perf);
    end
    % 
    if strcmp(obj.level,'subj') || ~plotError
        error = zeros(size(perf));
    else
        error = perfStruct.([perfIdx,'_',cvEval,'_err']);
        if ~isvector(error)
            error = diag(error);
        end
    end
    % Adding statistics if appropriate
    if addStats
        idx = ~cellfun(@isempty,regexp(statFields,...
                                       ['h_',fieldName],'once'));
        if any(idx)
            actH = stats{idx};
            clusterIdx = mvpares.findClusters(actH);
            det.addRect = arrayfun(@(x) time(x),clusterIdx);
        else
            warning('mvpares:showGenPerf:reqestedDataNotPresent',...
                    ['Couldn''t find statistics for this performance ',...
                     'measure in the dataset, plotting without statistics.']);
            det.addRect = [];
        end
    end
    % Preparing figure plotting
    if isfield(in,'f')
        det.hFig = in.f;
    end
    if isfield(in,'parentAxes')
        det.hAxes = in.parentAxes;
    end
    det.dataUnit = perfIdxStr;
    if isempty(lineProp)
        det.lineProp = {'Color','k','LineWidth',1.5,'LineStyle','-'};
    else
        det.lineProp = lineProp;
    end
    det.plotError = plotError;
    if ~isfield(in,'title')
        det.title = sprintf('Generalization performance over training time\n(train: %s, gen: %s)',...
            trCond,genCond);
    else
        det.title = in.title;
    end
    det.xLim = [sTime,eTime];
    hFig = mvpares.plotTimeSeriesData(time,perf,error,det);
    
else
    % Plotting decoding accuracy vs time.
    if strcmp(genTime,'tr')
        hFig = [];
        warning('mvpares:showGenPerf:reqestedDataNotPresent',...
                'The dataset is not generalized time x time. Returning');
        return;
    end
    
    % Adding statistics if appropriate
    if addStats
        idx = ~cellfun(@isempty,regexp(statFields,...
                                       ['h_',fieldName],'once'));
        if any(idx)
            actH = stats{idx};
            det.contour = actH;
        else
            warning('mvpares:showGenPerf:reqestedDataNotPresent',...
                    ['Couldn''t find statistics for this performance ',...
                     'measure in the dataset, plotting without statistics.']);
            det.contour = [];
        end
    end
    
    % Preparing figure plotting
    if isfield(in,'f')
        det.hFig = in.f;
    end
    if isfield(in,'parentAxes')
        det.hAxes = in.parentAxes;
    end
    if ~isfield(in,'title')
        det.title = sprintf('Generalization performance across time\n(train: %s, gen: %s)',...
            trCond,genCond);
    else
        det.title = in.title;
    end
    det.addCbar = true;
    if ismember(perfIdx,{'r','rf','b'})
        det.cLimSymm = true;
    end
    det.dataUnit = perfIdxStr;
    det.xLim = [sTime,eTime];
    mvpares.plotTimeByTimeMatrix(time,perf,det)
    
    % cLim = max([max(genPerf),abs(min(genPerf))]);
    % % Choosing the increment value between the displayed axis ticks
    % if tr_timePoints(end)-tr_timePoints(1) < 0.8
    %     incr = 0.1;
    % else
    %     incr = 0.2;
    % end
    % % Tick labels converted to milliseconds
    % xTickLab = round((tr_timePoints(1):incr:tr_timePoints(end))*1000);
    % yTickLab = round((gen_timePoints(1):incr:gen_timePoints(end))*1000);
    % % Positions of labels
    % xTickPos = find(ismember(tr_timePoints*1000,xTickLab));
    % yTickPos = find(ismember(gen_timePoints*1000,yTickLab));
    % % Plotting
    % imagesc(genPerf,'Parent',parentAxes);
    % colormap('jet');
    % c = colorbar;
    % set(parentAxes,'YTick',xTickPos,'YTickLabel',...
    %                cellfun(@num2str,num2cell(xTickLab),'UniformOutput',false));
    % set(parentAxes,'YDir','normal');
    % set(parentAxes,'XTick',yTickPos,'XTickLabel',...
    %                cellfun(@num2str,num2cell(yTickLab),'UniformOutput',false));
    % xlabel(parentAxes,'Generalization time (ms)');
    % ylabel(parentAxes,'Training time (ms)');
    % ylabel(c,perfIdxStr);
    % caxis([-cLim,cLim]);
    % % Figure title
    % if ~isfield(in,'title')
    %     title(parentAxes,...
    %           sprintf('Generalization performance across time\n(train: %s, gen: %s)',...
    %                   trCond,genCond));
    % else
    %     title(parentAxes,in.title);
    % end

    
    
    % genPerf = genPerf(genTimePointIdx);
    % yTickLab = round(gen_timePoints*1000);
    % plot(parentAxes,yTickLab,genPerf); hold on;
    % xlabel(parentAxes,'Time (ms)');
    % ylabel(parentAxes,perfIdxStr);
    
    % if ~isfield(in,'title')
    %     title(parentAxes,sprintf('Generalization performance over training time\n(train: %s, gen: %s)',...
    %                              trCond,genCond));
    % else
    %     title(parentAxes,in.title);
    % end

end
end