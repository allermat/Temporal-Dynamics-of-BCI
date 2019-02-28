function hFig = showAVmodelEstimates(obj,effect,varargin)
% Method for plotting AV model estimates
%
% USAGE:
%   hFig = showAVmodelEstimates(obj,effect,varargin)
%   hFig = showAVmodelEstimates(obj,effect,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%       effect (string): indicating which effects to plot. Possible values:
%           'int3way','int2way','main' respectively for 3-way interaction,
%           2-way interactions, main effects. Default: 'main'
%   'Name'-Value arguments:
%       plotError (logical): wheter to add error bars (takes effect only if
%           time courses are plotted (e.g. genTime is 'tr'). Default: true
%       smooth (logical): whether to smooth the time course of AV
%           estimates. Default: false
%       cust (struct): various settings to customize the plot's
%           appearence
%           Possible fields: supTitle, betaLims
% OUTPUT:
%   hFig (scalar/vector): figure handles for plotted figures 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validEffects = {'int3way','int2way','main'};
addRequired(p,'obj');
addRequired(p,'effect',@(x)any(validatestring(x,validEffects)));
addParameter(p,'plotError',true,@(x)validateattributes(x,{'logical'},...
                                                  {'scalar'}));
addParameter(p,'smooth',false,@(x)validateattributes(x,{'logical'},...
                                                  {'scalar'}));
addParameter(p,'cust',[],@isstruct);

parse(p,obj,effect,varargin{:});

obj = p.Results.obj;
effect = p.Results.effect;
plotError = p.Results.plotError;
smooth = p.Results.smooth;
cust = p.Results.cust;

estimates = obj.getAVmodelEstimates('smooth',smooth);
if isempty(estimates)
    warning('mvpares:showAVmodelEstimates:reqestedDataNotPresent',...
            'Couldn''t find AV model estimates in the dataset, returning.');
    return;
end
estimateFnames = fieldnames(estimates);

hFig = [];

time = obj.getTrTimePoints;
sTime = obj.plotTimeWin(1);
eTime = obj.plotTimeWin(2);
if isfield(cust,'betaLims') && ~isempty(cust.betaLims)
    betaLims = cust.betaLims;
else
    betaLims = [-0.2,0.6];
end

switch effect
  case 'int3way'
    hActFig = figure();
    hFig = [hFig,hActFig];
    det.hFig = hActFig;
    det.dataUnit = 'beta';
    condsOfInterest = {'A_r_d_a','A_r_d_v','A_R_d_a','A_R_d_v', ...
                       'A_r_D_a','A_r_D_v','A_R_D_a','A_R_D_v'; ...
                       'V_r_d_a','V_r_d_v','V_R_d_a','V_R_d_v', ...
                       'V_r_D_a','V_r_D_v','V_R_D_a','V_R_D_v'};
    % Checking which conditions are present in the dataset
    isCondPresent = true(size(condsOfInterest));
    for i = 1:numel(condsOfInterest)
        isCondPresent(i) = any(~cellfun(@isempty, ...
                                        regexp(estimateFnames,['^',condsOfInterest{i},'_beta'])));
    end
    % The length of this vector will be size(isCondPresent,1)
    if any(~isCondPresent(:))
        warning('mvpares:showAVmodelEstimates:reqestedDataNotPresent',...
                ['Couldn''t find AV model estimates required for ' ...
                 '3-way interactions. Returning']);
        return;
    end
    % condsOfInterest = condsOfInterest(isCondPresent);
    % The array needs to be rashaped to an N x 2 matrix
    % condsOfInterest = reshape(condsOfInterest,[],2);
    % Defining various variables related to the effects
    titles = {'Beta A','Beta V'};
    legends = {...
        'VR-, D-, A','VR-, D-, V','VR+, D-, A','VR+, D-, V', ...
        'VR-, D+, A','VR-, D+, V','VR+, D+, A','VR+, D+, V', ...
        'Location','SouthEastOutside';...
        'VR-, D-, A','VR-, D-, V','VR+, D-, A','VR+, D-, V', ...
        'VR-, D+, A','VR-, D+, V','VR+, D+, A','VR+, D+, V', ...
        'Location','SouthEastOutside'};
    lineProps = cat(1,...
                    {'Color','k','LineWidth',1.5,'LineStyle','-'},...
                    {'Color','k','LineWidth',1.5,'LineStyle','--'},...
                    {'Color','k','LineWidth',1.5,'LineStyle',':'},...
                    {'Color','k','LineWidth',1.5,'LineStyle','-.'},...
                    {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-'},...
                    {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','--'},...
                    {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle',':'},...
                    {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-.'});
    % titles = titles(isEffectPresent);
    % effectNames = effectNames(isEffectPresent);
    % legends = legends(isEffectPresent,:);
    % lineProps = lineProps(isEffectPresent);
    for i = 1:numel(titles)                
        actConds = condsOfInterest(i,:);
        actBetas = strcat(actConds,'_beta');
        actCis = strcat(actConds,'_ci');
        data = cat(2, ...
                   estimates.(actBetas{1}),estimates.(actBetas{2}), ...
                   estimates.(actBetas{3}),estimates.(actBetas{4}), ...
                   estimates.(actBetas{5}),estimates.(actBetas{6}), ...
                   estimates.(actBetas{7}),estimates.(actBetas{8}));
        error = cat(2, ...
                    estimates.(actCis{1}),estimates.(actCis{2}), ...
                    estimates.(actCis{3}),estimates.(actCis{4}), ...
                    estimates.(actCis{5}),estimates.(actCis{6}), ...
                    estimates.(actCis{7}),estimates.(actCis{8}));
        det.lineProp = lineProps;
        det.hAxes = subplot(1,numel(titles),i);
        det.legend = legends(i,:);
        det.plotError = plotError;
        det.title = titles{i};
        det.yLim = betaLims;
        det.xLim = [sTime,eTime];
        mvpares.plotTimeSeriesData(time,data,error,det);
    end
    % Setting figure size
    set(hActFig,'Units','normalized','Position',[1/6,1/4,2/3,1/2],...
                'PaperPositionMode','auto');
    if isfield(cust,'supTitle')
        suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
    else
        suplabel('AV model estimates, 3-way interactions','t',[0.08,0.08,0.84,0.88]);
    end    
    
  case 'int2way'
    hActFig = figure();
    hFig = [hFig,hActFig];
    det.hFig = hActFig;
    det.dataUnit = 'beta';
    condsOfInterest = {...
        'A_r_d','A_r_D','A_R_d','A_R_D'; ...
        'V_r_d','V_r_D','V_R_d','V_R_D'; ...
        'A_r_a','A_r_v','A_R_a','A_R_v'; ...
        'V_r_a','V_r_v','V_R_a','V_R_v'; ...
        'A_D_a','A_d_v','A_D_a','A_D_v'; ...
        'V_d_a','V_d_v','V_D_a','V_D_v'};
    isCondPresent = true(size(condsOfInterest));
    for i = 1:numel(condsOfInterest)
        isCondPresent(i) = any(~cellfun(@isempty, ...
                                        regexp(estimateFnames,['^',condsOfInterest{i},'_beta'])));
    end
    % The length of this vector will be size(isCondPresent,1)
    isEffectPresent = all(isCondPresent,2);
    % condsOfInterest = condsOfInterest(isCondPresent);
    % The array needs to be rashaped to an N x 4 matrix
    % condsOfInterest = reshape(condsOfInterest,[],4);
    titles = {'VR X Disparity, beta A', 'VR X Disparity, beta V', ...
              'VR X Task, beta A','VR X Task, beta V', ...
              'Disparity X Task, beta A','Disparity X Task, beta V'};
    legends = {...
        'VR-, D-','VR-, D+','VR+, D-','VR+, D+','Location','SouthEastOutside';...
        'VR-, D-','VR-, D+','VR+, D-','VR+, D+','Location','SouthEastOutside';...
        'VR-, A','VR-, V','VR+, A','VR+, V','Location','SouthEastOutside';...
        'VR-, A','VR-, V','VR+, A','VR+, V','Location','SouthEastOutside';...
        'D-, A','D-, V','D+, A','D+, V','Location','SouthEastOutside'; ...
        'D-, A','D-, V','D+, A','D+, V','Location','SouthEastOutside'};
    % titles = titles(isEffectPresent);
    % legends = legends(isEffectPresent,:);
    for i = 1:numel(titles)
        if ~isEffectPresent(i), continue; end
        actConds = condsOfInterest(i,:);
        actWavs = strcat(actConds,'_beta');
        actCis = strcat(actConds,'_ci');
        data = cat(2,...
                   estimates.(actWavs{1}),estimates.(actWavs{2}),...
                   estimates.(actWavs{3}),estimates.(actWavs{4}));
        error = cat(2,...
                    estimates.(actCis{1}),estimates.(actCis{2}),...
                    estimates.(actCis{3}),estimates.(actCis{4}));
        det.lineProp = cat(1,...
                           {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-'},...
                           {'Color','k','LineWidth',1.5,'LineStyle','-'},...
                           {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','--'},...
                           {'Color','k','LineWidth',1.5,'LineStyle','--'});
        det.hAxes = subplot(numel(titles)/2,2,i);
        det.legend = legends(i,:);
        det.plotError = plotError;
        det.title = titles{i};
        det.yLim = betaLims;
        det.xLim = [sTime,eTime];
        det.addLine = {[[sTime eTime]' [sTime eTime]'],[[90 90]' [0 0]'],...
                       'Color',[0.5 0.5 0.5]};
        mvpares.plotTimeSeriesData(time,data,error,det);
    end
    % Setting figure size
    set(hActFig,'Units','normalized','Position',[1/6,0,2/3,1],...
                'PaperPositionMode','auto');
    if isfield(cust,'supTitle')
        suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
    else
        suplabel('AV model estimates, 2-way interactions','t',[0.08,0.08,0.84,0.88]);
    end
  case 'main'
    hActFig = figure();
    hFig = [hFig,hActFig];
    det.hFig = hActFig;
    det.dataUnit = 'beta';
    condsOfInterest = {'A_r','A_R','V_r','V_R'; ...
                       'A_a','A_v','V_a','V_v'; ...
                       'A_d','A_D','V_d','V_D'};
    % Checking which conditions are present in the dataset
    isCondPresent = true(size(condsOfInterest));
    for i = 1:numel(condsOfInterest)
        isCondPresent(i) = any(~cellfun(@isempty, ...
                                        regexp(estimateFnames,['^',condsOfInterest{i},'_beta'])));
    end
    % The length of this vector will be size(isCondPresent,1)
    isEffectPresent = all(isCondPresent,2);
    % condsOfInterest = condsOfInterest(isCondPresent);
    % The array needs to be rashaped to an N x 2 matrix
    % condsOfInterest = reshape(condsOfInterest,[],2);
    % Defining various variables related to the effects
    titles = {'Visual reliability','Task','Disparity'};
    effectNames = {'VR','Task','Disp'};
    legends = {...
        'bA, VR-','bA, VR+','bV, VR-','bV, VR+','Location','SouthEastOutside';...
        'bA, A','bA, V','bV, A','bV, V','Location','SouthEastOutside';...
        'bA, D-','bA, D+','bV, D-','bV, D+','Location','SouthEastOutside'};
    lineProps = {...
        cat(1,...
            {'Color','k','LineWidth',1.5,'LineStyle','-'},...
            {'Color','k','LineWidth',1.5,'LineStyle','--'},...
            {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-'},...
            {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','--'}),...
        cat(1,...
            {'Color','k','LineWidth',1.5,'LineStyle','-'},...
            {'Color','k','LineWidth',1.5,'LineStyle','--'},...
            {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-'},...
            {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','--'}),...
        cat(1,...
            {'Color','k','LineWidth',1.5,'LineStyle','-'},...
            {'Color','k','LineWidth',1.5,'LineStyle','--'},...
            {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','-'},...
            {'Color',[0.5 0.5 0.5],'LineWidth',1.5,'LineStyle','--'})};
    % titles = titles(isEffectPresent);
    % effectNames = effectNames(isEffectPresent);
    % legends = legends(isEffectPresent,:);
    % lineProps = lineProps(isEffectPresent);
    for i = 1:numel(effectNames)                
        % Simply skip those effects which are not present
        if ~isEffectPresent(i), continue; end
        actConds = condsOfInterest(i,:);
        actBetas = strcat(actConds,'_beta');
        actCis = strcat(actConds,'_ci');
        data = cat(2,estimates.(actBetas{1}),estimates.(actBetas{2}), ...
                   estimates.(actBetas{3}),estimates.(actBetas{4}));
        error = cat(2,estimates.(actCis{1}),estimates.(actCis{2}), ...
                    estimates.(actCis{3}),estimates.(actCis{4}));
        det.lineProp = lineProps{i};
        det.hAxes = subplot(numel(effectNames),1,i);
        det.legend = legends(i,:);
        det.plotError = plotError;
        det.title = titles{i};
        det.yLim = betaLims;
        det.xLim = [sTime,eTime];
        mvpares.plotTimeSeriesData(time,data,error,det);
    end
    % Setting figure size
    set(hActFig,'Units','normalized','Position',[1/3,0,1/3,1],...
                'PaperPositionMode','auto');
    if isfield(cust,'supTitle')
        suplabel(cust.supTitle,'t',[0.08,0.08,0.84,0.88]);
    else
        suplabel('AV model estimates, main effects','t',[0.08,0.08,0.84,0.88]);
    end    
    
end


end