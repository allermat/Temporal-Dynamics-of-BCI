function hFig = showPredLabels(obj,varargin)
% Method for plotting predicted labels vs. true labels
%
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       poolFolds (logical): whether to plot the predicted labels pooled
%           over CV folds, default: true
%       zScore (logical): whether to Z-score the true labels prior to
%           plotting, default: false
% OUTPUT:
%   hFig (scalar/vector): figure handles for plotted figures
%
% Copyright(C) 2016, Mate Aller

% Parsing input
p = inputParser;
addRequired(p,'obj');
addParameter(p,'poolFolds',true,@(x)validateattributes(x,...
    {'logical'},{'scalar','nonempty'}));
addParameter(p,'zScore',false,@(x)validateattributes(x,...
    {'logical'},{'scalar','nonempty'}));
parse(p,obj,varargin{:});
obj = p.Results.obj;
poolFolds = p.Results.poolFolds;
zScore = p.Results.zScore;

% Each figure contains a number of subplots
nRowsFig = 5;
nColsFig = 7;
nSubplots = nRowsFig*nColsFig;
% The time points which will be represented are chosen
trTimePoints = obj.getTrTimePoints;
n = ceil(numel(trTimePoints)/nSubplots);
timePointIdx = 1:n:numel(trTimePoints);
timePointsToPlot = trTimePoints(timePointIdx);

% Plotting actual locations vs decoded locations for each timepoint.
if poolFolds
    % Pooling over models
    predLabels = obj.getPredLabels('trTimePoints',timePointsToPlot);
    trueLabels = obj.getTrueLabels;
    if any([isempty(predLabels),isempty(trueLabels)])
        hFig = [];
        warning('mvpares:showPredLabels:missingLabels',...
            'True or predicted labels are missing');
        return;
    end
    % Z-scoring, if necessary
    if zScore
        m = nanmean(predLabels(:));
        s = nanstd(predLabels(:));
        predLabels = (predLabels-m)/s;
    end
    % Finding the max and min of the data
    y_lim = prctile(predLabels(:),[0.1,99.9]);
    % Plotting figure
    hFig  = figure();
    for i = 1:numel(timePointsToPlot)
        h = subplot(nRowsFig,nColsFig,i);
        if i/nColsFig > nRowsFig-1
            xLabelStr = 'Actual location';
        else
            xLabelStr = '';
        end
        if mod(i,nColsFig) == 1
            yLabelStr = 'Decoded location';
        else
            yLabelStr = '';
        end
        titleStr = sprintf('%d ms',timePointsToPlot(i)*1000);
        obj.plotLabelDistribution(h,predLabels(:,i),trueLabels,...
            xLabelStr,yLabelStr,titleStr,y_lim)
    end
    suplabel(sprintf('Actual locations vs. decoded locations (%s).',...
        obj.info.tr_cond),'t',[0.08 0.08 0.88 0.88]);
    % Setting figure size
    set(hFig,'Units','normalized','Position',[0,0,1,1],...
        'PaperPositionMode','auto');
else
    % There will be one figure for each fold
    nPlots = obj.getNcvFolds;
    % Vector for figure handles
    hFig = NaN(1,nPlots);
    for j = 1:nPlots
        predLabels = obj.getPredLabels('cvFold',j,'trTimePoints',timePointsToPlot);
        trueLabels = obj.getTrueLabels('cvFold',j);
        if any([isempty(predLabels),isempty(trueLabels)])
            hFig = [];
            warning('mvpares:showPredLabels:missingLabels',...
                'True or predicted labels are missing');
            return;
        end
        % Z-scoring, if necessary
        if zScore
            m = nanmean(predLabels(:));
            s = nanstd(predLabels(:));
            predLabels = (predLabels-m)/s;
        end
        % Finding the max and min of the data
        y_lim = prctile(predLabels(:),[0.1,99.9]);
        % Plotting figures
        hActFig = figure();
        hFig(j) = hActFig;
        for i = 1:numel(timePointsToPlot)
            h = subplot(nRowsFig,nColsFig,i);
            if i/nColsFig > nRowsFig-1
                xLabelStr = 'Actual location';
            else
                xLabelStr = '';
            end
            if mod(i,nColsFig) == 1
                yLabelStr = 'Decoded location';
            else
                yLabelStr = '';
            end
            titleStr = sprintf('%d ms',timePointsToPlot(i)*1000);
            obj.plotLabelDistribution(h,predLabels(:,j),trueLabels,...
                xLabelStr,yLabelStr,titleStr,y_lim);
        end
        suplabel(...
            sprintf('Actual locations vs. decoded locations (%s), fold %d.',...
            obj.info.tr_cond,j),'t',[0.08 0.08 0.88 0.88]);
        % Setting figure size
        set(hActFig,'Units','normalized','Position',[0,0,1,1],...
            'PaperPositionMode','auto');
    end
end

end