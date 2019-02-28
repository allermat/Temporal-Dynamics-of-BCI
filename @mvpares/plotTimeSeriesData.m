function [hFig,hAxes] = plotTimeSeriesData(time,data,varargin)
% Method for plotting data of timeseries format
%
% USAGE:
%   [hFig,hAxes] = plotTimeSeriesData(time,data)
%   [hFig,hAxes] = plotTimeSeriesData(time,data,error)
%   [hFig,hAxes] = plotTimeSeriesData(time,data,error,det)
% INPUT:
%   Required:
%       time (numeric vector): time values in seconds
%       data (numeric vector): data to be plotted.
%   Optional:
%       error (numeric vector): error values to be plotted with data. 
%       det (struct): Additional settings and input as a structure.
%           Possible fields: 
%               det.addLine
%               det.addRect
%               det.dataUnit
%               det.hFig
%               det.hAxes
%               det.legend
%               det.lineProp
%               det.plotError
%               det.title
%               det.transp
%               det.xLim
%               det.yLim
%               det.yLimSymm
%
% OUTPUT:
%   hFig (scalar): figure handle for plotted figure
%   hAxes (scalar): axes handle for plotted figure

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input, checking matlab
p = inputParser;
addRequired(p,'time',@(x) validateattributes(x,{'numeric'},{'vector',...
    'numel',size(data,1),'increasing'}));
addRequired(p,'data',@(x) validateattributes(x,{'numeric'},{'ndims',2,...
    'nonempty'}));
addOptional(p,'error',[],@(x) validateattributes(x,{'numeric'},{'size',size(data)}));
addOptional(p,'det',struct(),@(x) validateattributes(x,{'struct'},{'nonempty'}));
parse(p,time,data,varargin{:});
data = p.Results.data;
time = p.Results.time;
error = p.Results.error;
det = p.Results.det;

% Preparing figure for plotting
if isfield(det,'hFig')
    hFig = det.hFig;
else
    hFig = figure();
end
if isfield(det,'hAxes')
    hAxes = det.hAxes;
else
    hAxes = axes('parent',hFig);
end
if isfield(det,'lineProp'), lineProp = det.lineProp; else lineProp = {}; end
if isfield(det,'transp'), transp = det.transp; else transp = true; end
if isempty(error)
    plotError = false;
    if isfield(det,'plotError') && det.plotError
        warning('mvpares:plotTimeSeriesData:reqestedDataNotPresent',...
            'Error estimates were not provided, plotting without error bars. ');
    end
else
    if isfield(det,'plotError') && ~det.plotError
        plotError = false;
    else
        plotError = true;
    end
end
   
% Converting seconds to milliseconds
time = round(time*1000);

% Plotting
if size(data,2) == 1
    if plotError
        h = shadedErrorBar(time,data,error,lineProp,transp);
    else
        plot(hAxes,time,data,lineProp{:});
    end
else
    if plotError
        for i = 1:size(data,2)
            if isempty(lineProp), 
                actLineProp = {};
            else
                actLineProp = lineProp(i,:);
            end
            h(i) = shadedErrorBar(time,data(:,i),error(:,i),actLineProp,transp); %#ok
            if i == 1, hold on; end
        end
    else
        for i = 1:size(data,2)
            if isempty(lineProp)
                actLineProp = {};
            else
                actLineProp = lineProp(i,:);
            end
            plot(hAxes,time,data(:,i),actLineProp{:});
            if i == 1, hold on; end
        end
    end
end
% xlim
if isfield(det,'xLim')
    xlim(hAxes,det.xLim*1000);
else
    xlim(hAxes,[min(time),max(time)]);
end
% ylim
if isfield(det,'yLim'), ylim(hAxes,det.yLim); end
% ylimSymm
if isfield(det,'yLimSymm') && det.yLimSymm
    yl = ylim;
    yl = max(abs(yl));
    ylim([-yl,yl]);
end
% Add line(s) if necessary
if isfield(det,'addLine')
    X = det.addLine{1}*1000;
    Y = det.addLine{2};
    misc = det.addLine(3:end);
    line(X,Y,misc{:});
end
% Add rectangle(s) if necessary
if isfield(det,'addRect')
    rect = det.addRect*1000;
    yl = ylim;
    for i = 1:size(rect,1)
        hr = rectangle('Position',[rect(i,1),yl(1),abs(rect(i,2)-rect(i,1)),abs(yl(2)-yl(1))],...
            'Parent',hAxes,'LineStyle','none','FaceColor',[0.9 0.9 0.9]);
        uistack(hr,'bottom');
    end
end
% x axis label
xlabel(hAxes,'Training time (ms)');
% Legend
if isfield(det,'legend')
    if plotError
        legend([h.mainLine],det.legend{:});
    else
        legend(det.legend{:});
    end
end
% Unit of data
if isfield(det,'dataUnit'), ylabel(hAxes,det.dataUnit); end

% Figure title
if isfield(det,'title'), title(hAxes,det.title); end


end

