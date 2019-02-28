function [hFig,hAxes] = plotTimeByTimeMatrix(time,data,varargin)
% Method for plotting data of time x time matrix format
%
% USAGE:
%   [hFig,hAxes] = plotTimeByTimeMatrix(time,data)
%   [hFig,hAxes] = plotTimeByTimeMatrix(time,data,det)
% INPUT:
%   Required:
%       time (numeric vector): time values in seconds
%       data (numeric square matrix): data to be plotted.
%   Optional:
%       det (struct): Additional settings and input as a structure.
%           Possible fields: 
%               det.addLines
%               det.addCbar
%               det.cLim
%               det.cLimSymm
%               det.dataUnit
%               det.hFig
%               det.hAxes
%               det.lineCol
%               det.title
%               det.xLab
%               det.xLim
%               det.yLab
%
% OUTPUT:
%   hFig (scalar): figure handle for plotted figure
%   hAxes (scalar): axes handle for plotted figure
%
% Copyright(C) 2016, Mate Aller

% Parsing input
p = inputParser;
addRequired(p,'time',@(x) validateattributes(x,{'numeric'},{'vector',...
    'numel',size(data,1),'increasing'}));
addRequired(p,'data',@(x) validateattributes(x,{'numeric'},{'square',...
    'nonempty'}));
addOptional(p,'det',struct(),@(x) validateattributes(x,{'struct'},{'nonempty'}));
parse(p,time,data,varargin{:});
data = p.Results.data;
time = p.Results.time;
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
% Finding data and timepoints of interest
if isfield(det,'xLim')
    winOfIntIdx = time >= det.xLim(1) & time <= det.xLim(2);
    time = time(winOfIntIdx);
    data = data(winOfIntIdx,winOfIntIdx);
end
% Choosing the increment value between the displayed axis ticks
if time(end)-time(1) < 0.8
    incr = 0.1;
else
    incr = 0.2;
end
% Tick labels converted to milliseconds
xTickLab = round((time(1):incr:time(end))*1000);
yTickLab = xTickLab;
% Positions of labels
xTickPos = find(ismember(time*1000,xTickLab));
yTickPos = find(ismember(time*1000,yTickLab));

% Plotting
imagesc(data,'Parent',hAxes); hold on;

% Setting colormap
colormap('redblue');
if ~isfield(det,'addLines') || det.addLines
    if isfield(det,'lineCol'), lineCol = det.lineCol; else lineCol = 'k'; end
    if any(time == 0)
        line([find(time == 0),find(time == 0)],[0,numel(time)],'Color',lineCol);
        line([0,numel(time)],[find(time == 0),find(time == 0)],'Color',lineCol);
    end
    line([0,numel(time)],[0,numel(time)],'Color',lineCol);
end
if isfield(det,'addCbar') && det.addCbar
    hCaxis = colorbar;
end

set(hAxes,'YDir','normal');
set(hAxes,'YTick',yTickPos,'YTickLabel',yTickLab);
set(hAxes,'XTick',xTickPos,'XTickLabel',xTickLab);
if isfield(det,'xLab')
    xlabel(hAxes,det.xLab);
else
    xlabel(hAxes,'Generalization time (ms)');
end
if isfield(det,'yLab')
    ylabel(hAxes,det.yLab);
else
    ylabel(hAxes,'Training time (ms)');
end
% Unit of data
if isfield(det,'dataUnit') && isfield(det,'addCbar') && det.addCbar
    ylabel(hCaxis,det.dataUnit);
end

% Colorbar limits
if isfield(det,'cLim'), caxis(det.cLim); end

if isfield(det,'cLimSymm') && det.cLimSymm
    c = caxis;
    c = max(abs(c));
    caxis([-c,c]);
end

% Adding contour plot
if isfield(det,'contour') && ~isempty(det.contour)
    contour(det.contour,[1,1],'k');
end

% Figure title
if isfield(det,'title'), title(hAxes,det.title); end

% Making axes square
axis('square');

end

