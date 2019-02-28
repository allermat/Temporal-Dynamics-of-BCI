function [nRowSubplot,nColSubplot,subPlotIdx] = segmentFigArea(nRows,nCols,spacing)
% Method for segmenting figure area for multiple plots
% 
% USAGE: 
%   [nRowSubplot,nColSubplot,subPlotIdx] = segmentFigArea(nRows,nCols,spacing) 
% INPUT:
%   nRows (scalar): number of rows of subplots on figure
%   nCols (scalar): number of columns of subplots on figure
%   spacing (scalar): specifies the width of the space between subplots
%       relative to the width/height of the subplots. Varies between 0
%       and 1
% OUTPUT:
%   nRowSubplot (scalar): number of rows of subplots to be passed to the
%       subplot function to achieve the desired segmentation.
%   nColSubplot (scalar): number of columns of subplots to be passed to the
%       subplot function to achieve the desired segmentation.
%   subPlotIdx (cell): indices of subplots for plotting each subplot. Each
%       cell contains one or more indices which define the area on the
%       figure where the given subplot is to be drawn. 
% 
% Copyright(C) 2016, Mate Aller

% Parsing input
p = inputParser;
addRequired(p,'nRows',@(x) validateattributes(x,{'numeric'},{'scalar',...
    'integer','positive'}));
addRequired(p,'nCols',@(x) validateattributes(x,{'numeric'},{'scalar',...
    'integer','positive'}));
addRequired(p,'spacing',@(x) validateattributes(x,{'numeric'},{'scalar',...
    'nonnegative','<=',1}));
parse(p,nRows,nCols,spacing);
nRows = p.Results.nRows;
nCols = p.Results.nCols;
spacing = p.Results.spacing;

if spacing > 0
    nRowSubplot = round(nRows/spacing)+nRows+1;
    nColSubplot = round(nCols/spacing)+nCols+1;
else
    nRowSubplot = nRows;
    nColSubplot = nCols;
    subPlotIdx = num2cell((1:nRows*nCols)');
    return;
end

gridRowIdx = cell(nRows,1);
gridColIdx = cell(1,nCols);
nUnit = round(1/spacing)+2;
for i = 1:nRows
    temp = (2:(nUnit-1))+((i-1)*(nUnit-1));
    temp = repmat(temp,round(1/spacing),1);
    gridRowIdx{i} = temp(:);
end
for i = 1:nCols
    temp = (2:(nUnit-1))'+((i-1)*(nUnit-1));
    temp = repmat(temp,1,round(1/spacing));
    gridColIdx{i} = temp(:);
end

gridRowIdx = repmat(gridRowIdx,1,nCols);
gridColIdx = repmat(gridColIdx,nRows,1);
gridSize = repmat({[nColSubplot,nRowSubplot]},nRows,nCols);

subPlotIdx = cellfun(@sub2ind,gridSize,gridColIdx,gridRowIdx,'UniformOutput',false);
subPlotIdx = subPlotIdx';
subPlotIdx = subPlotIdx(:);

end

