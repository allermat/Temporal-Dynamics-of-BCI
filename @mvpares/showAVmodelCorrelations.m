function hFig = showAVmodelCorrelations(obj,var,varargin)
% Method for plotting AV model correlations
%
% USAGE:
%   hFig = showAVmodelCorrelations(obj)
%   hFig = showAVmodelCorrelations(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%       var (string): the variable whos AV weights were correlated with 
%           the AV weights of the mvpares object. Possible values:
%           'acrossTime', 'behav', 'fmri'
%   'Name'-Value arguments:
%       addStats (logical): wheter to add markers indicating significant
%           effects. Default: false
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       smooth (logical): whether to smooth the time courses 
% OUTPUT:
%   hFig (scalar/vector): figure handles for plotted figures 
%
% Copyright(C) 2016, Mate Aller

% Parsing input
p = inputParser;
validVars = {'acrossTime','behav','fmri'};
validGenTimes = {'tr','tr_x_tr'};
addRequired(p,'obj');
addRequired(p,'var',@(x) ismember(x,validVars));
addParameter(p,'addStats',true,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'genTime','tr',@(x) any(validatestring(x,validGenTimes)));
addParameter(p,'smooth',false,@(x) validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'fisherTransform',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
parse(p,obj,var,varargin{:});
obj = p.Results.obj;
addStats = p.Results.addStats;
genTime = p.Results.genTime;
var = p.Results.var;
smooth = p.Results.smooth;
fisherTransform = p.Results.fisherTransform;

corrCoeffs = obj.getAVmodelCorrelations('all','genTime',genTime,'smooth',smooth,...
                                        'fisherTransform',fisherTransform);

if isempty(corrCoeffs)
    hFig = [];
    warning('mvpares:showAVmodelEstimates:reqestedDataNotPresent',...
        'Couldn''t find AV model correlations in the dataset, returning.');
    return;
end
if addStats
    stats = obj.getStats('avModelCorr','genTime',genTime,'smooth',smooth);
    if isempty(stats)
        warning('mvpares:showAVmodelWeights:reqestedDataNotPresent',...
            ['Couldn''t find statistics for AV model correlations ',...
             'in the dataset, plotting without statistics.']);
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
tWin = obj.plotTimeWin;

switch var
    case 'acrossTime'
        data = corrCoeffs.cc_acrossTime;
        det.title = 'Across time correlation of neural AV weights';
        if fisherTransform
            det.dataUnit = 'R (Fisher transformed)';
            det.cLimSymm = true;
        else
            det.dataUnit = 'R';
            det.cLim = [-1,1];
        end
        det.xLim = tWin;
        
        det.addCbar = true;
        det.xLab = 'time (ms)';
        det.yLab = 'time (ms)';
        if addStats
            idx = ~cellfun(@isempty,regexp(statFields,...
                                           'h_.*_acrossTime','once'));
            if any(idx)
                det.contour = stats{idx};
            else
                det.addRect = [];
            end
        end
        hFig = mvpares.plotTimeByTimeMatrix(time,data,det);

    case 'behav'
        data = corrCoeffs.cc_behav;
        error = corrCoeffs.ci_behav;
        if addStats
            if all(cellfun(@isempty,regexp(statFields,'h_.*_behav','once')))
                warning('mvpares:showAVmodelCorrelations:reqestedDataNotPresent',...
                    ['Couldn''t find statistics for AV model correlations for ',...
                    'this particular variable in the dataset, plotting ',...
                    'without statistics.']);
                addStats = false;
            end
        end
        det.title = 'Correlation between neural and behavioural AV weights';
        if fisherTransform
            det.dataUnit = 'R (Fisher transformed)';
        else
            det.dataUnit = 'R';
        end
        det.xLim = tWin;
        switch genTime
            case 'tr'
                det.lineProp = {'Color','k','LineWidth',1.5,'LineStyle','-'};
                if ~fisherTransform
                    det.yLim = [min(data),1];
                end
                if addStats
                    idx = ~cellfun(@isempty,regexp(statFields,...
                                                  'h_.*_behav','once'));
                    if any(idx)
                        actH = stats{idx};
                        clusterIdx = mvpares.findClusters(actH);
                        det.addRect = arrayfun(@(x) time(x),clusterIdx);
                    else
                        det.addRect = [];
                    end
                end
                hFig = mvpares.plotTimeSeriesData(time,data,error,det);
            case 'tr_x_tr'
                if fisherTransform
                    det.cLimSymm = true;
                else
                    det.cLim = [-1,1];
                end
                det.addCbar = true;
                hFig = mvpares.plotTimeByTimeMatrix(time,data,det);
        end
    case 'fmri'
        fields = fieldnames(corrCoeffs);
        fmriCorrNames = fields(~cellfun(@isempty,regexp(fields,'^cc_fmri','once')));
        % Reordering according to the hierarchy
        fmriCorrNames = fmriCorrNames([4:7,2:3,1,8]);
        temp = regexp(fmriCorrNames,'^cc_fmri_(.*)','tokens');
        rois = [temp{:}]';
        hFig = figure();
        det.hFig = hFig;
        nColFig = 4;
        nRowFig = ceil(numel(fmriCorrNames)/nColFig);
        switch genTime
            case 'tr'
                for i = 1:numel(fmriCorrNames)
                    data = corrCoeffs.(fmriCorrNames{i});
                    error = corrCoeffs.(strrep(fmriCorrNames{i},'cc_','ci_'));
                    det.hAxes = subplot(nRowFig,nColFig,i);
                    det.title = strrep(rois{i},'_','\_');
                    if fisherTransform
                        det.dataUnit = 'R (Fisher transformed)';
                        det.yLim = [min(data),max(data)];
                    else
                        det.dataUnit = 'R';
                        det.yLim = [-0.2,1];
                    end
                    det.lineProp = {'Color','k','LineWidth',1.5,'LineStyle','-'};
                    
                    det.xLim = tWin;
                    mvpares.plotTimeSeriesData(time,data,error,det);
                end
                set(hFig,'Units','normalized','Position',[1/8,1/4,3/4,0.6],...
                    'PaperPositionMode','auto');
                suplabel('Correlation between EEG and fMRI AV weights','t');
            case 'tr_x_tr'
                for i = 1:numel(fmriCorrNames)
                    data = corrCoeffs.(fmriCorrNames{i});
                    det.hAxes = subplot(nRowFig,nColFig,i);
                    det.title = strrep(rois{i},'_','\_');
                    if fisherTransform
                        det.cLimSymm = true;
                    else
                        det.cLim = [-1,1];
                    end
                    det.xLim = tWin;
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
                if fisherTransform
                    ylabel(hCbar,'R (Fisher transformed)');
                else
                    ylabel(hCbar,'R');
                end
                
                suplabel('Correlation between EEG and fMRI AV weights','t');
        end
end

end