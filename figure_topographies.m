function figure_topographies(subID,time,conds,varargin)
% Plots a series of topographies for four locations at specified time points
% 
%% Parsing input, checking matlab
p = inputParser;

validConds = {'A','V_rh','V_rl','AV_rh_a','AV_rl_a','AV_rh_v','AV_rl_v'};
validOrientations = {'horiz','vert'};

addRequired(p,'subID',@ischar);
addRequired(p,'time',@isnumeric);
addRequired(p,'conds',@(x)all(ismember(x,validConds)));
% addOptional(p,'weights',[],@isnumeric);
% addOptional(p,'weight_times',[],@isnumeric);
addParameter(p,'orientation','horiz',@(x)any(validatestring(x,validOrientations)))
addParameter(p,'averageConds',[],@(x) validateattributes(x,{'numeric'},...
    {'vector','>=',2,'<=',numel(conds)}))

parse(p,subID,time,conds,varargin{:});

subID = p.Results.subID;
time = p.Results.time;
conds = p.Results.conds;
% weights = p.Results.weights;
% weight_times = p.Results.weight_times;
orientation = p.Results.orientation;
averageConds = p.Results.averageConds;

%% Loading data files
if strcmp(subID,'group')
    dataFolder = fullfile(DEC_2_setupdir('final','anal_eeg_group'),'ERP');
else
    dataFolder = DEC_2_setupdir('final','anal_eeg_sub_erp',subID);
end
saveDf = cd(dataFolder);
fileList = cellstr(ls('*.mat'));
cd(saveDf);

for i = 1:numel(conds)
    % Defining matchstrings for conditions
    if strcmp(conds{i},'A')
        matchStr = '_A_([-0-9]*)';
    elseif strcmp(conds{i},'V_rh')
        matchStr = '_V_([-0-9]*)_rh';
    elseif strcmp(conds{i},'V_rl')
        matchStr = '_V_([-0-9]*)_rl';
    elseif strcmp(conds{i},'AV_rh_a')
        matchStr = '_AV_([-0-9]*)_rh_a';
    elseif strcmp(conds{i},'AV_rl_a')
        matchStr = '_AV_([-0-9]*)_rl_a';
    elseif strcmp(conds{i},'AV_rh_v')
        matchStr = '_AV_([-0-9]*)_rh_v';
    elseif strcmp(conds{i},'AV_rl_v')
        matchStr = '_AV_([-0-9]*)_rl_v';
    end
    
    found = regexp(fileList,matchStr,'match','once');
    found = ~strcmp(found,'');
    if all(~found)
        error('No files match the specified condition!');
    end
    actSourceFileNames = fileList(found);
    loc = regexp(fileList,matchStr,'tokens');
    loc = loc(found);
    loc = [loc{:}]';
    loc = vertcat(loc{:});
    loc = cell2mat(cellfun(@str2num,loc,'UniformOutput',false));
    [loc,locOrder] = sort(loc);
    % Putting the source files in order
    actSourceFileNames = actSourceFileNames(locOrder);
    sourceFileNames(i,:) = actSourceFileNames;
end

% Pre allocating arrays for collecting data
ftDataStructures = cell(size(sourceFileNames));

%% 
% Loading the eeg data
allData = [];
for i = 1:numel(sourceFileNames)
    temp = load(fullfile(dataFolder,sourceFileNames{i}));
    tempFieldName = fieldnames(temp);
    ftDataStructures{i} = temp.(tempFieldName{1});
    allData = cat(3,allData,ftDataStructures{i}.avg);
end

if ~isempty(averageConds)
    temp = {};
    for i = 1:size(ftDataStructures,2)
        temp{1,i} = ft_appenddata([],ftDataStructures{averageConds,i});
        temp{1,i} = ft_timelockanalysis([],temp{1,i});
        temp{1,i}.avg = temp{1,i}.avg.*2;
    end
    ftDataStructures(averageConds,:) = [];
    ftDataStructures = cat(1,ftDataStructures,temp);
    temp = strcat(conds(averageConds),'+');
    temp = strcat(temp{:});
    conds(averageConds) = [];
    conds = cat(2,conds,temp(1:end-1));
end

% Checking the maximum and minimum of the used data
timeIdx = false(size(ftDataStructures{1}.time));
for i = 1:numel(time)
    idx = find(abs(ftDataStructures{1}.time-time(i)) <= eps('double'));
    if ~isempty(idx), timeIdx(idx) = true; end    
end
allData = allData(:,timeIdx);
minV = min(allData(:));
maxV = max(allData(:));

nRowSubplots = size(ftDataStructures,1);
nColSubplots = size(ftDataStructures,2);
nSubplots = numel(ftDataStructures);

% The subplot function fills the rows first, hence we transpose the
% data matrix
ftDataStructures = ftDataStructures';

for iTimePoint = 1:numel(time)
    
    in.f = figure('Units','Normalized','Position',[0.2,0.125,0.6,0.75]);
    in.noButtons = true;
    in.cbar = 0;
    
%     if ~isempty(weights)
%         
%         if isempty(weight_times)
%             error('weight_times must be specified if weights are specified');
%         end
%         wSample = find(weight_times == time(iTimePoint));
%         if isempty(wSample)
%             error('The specified time point is not present in the weights');
%         end
%         nSubplots = nSubplots+1;
%         w = squeeze(mean(weights,2));
%         minW = min(w(:));
%         maxW = max(w(:));
%         
%     end
    
    % h = cell(1,numel(conds));
    
    for iPlot = 1:nSubplots
        
        % h{iPlot} = subplot(nRowSubplots,nColSubplots,iPlot);
        switch orientation
            case 'horiz'
                h{iPlot} = subplot(nRowSubplots,nColSubplots,iPlot);
            case 'vert'
                h{iPlot} = subplot(nSubplots,1,iPlot);
        end
        in.ParentAxes = h{iPlot};
        
        
%         if ~isempty(weights) && iPlot == nSubplots
%             data = w(:,wSample);
%             titleStr = 'Feature weights';
%         else
%             data = ftDataStructures{iPlot};
%             titleStr = num2str(loc(iPlot));
%         end
        data = ftDataStructures{iPlot};
        titleStr = num2str(loc(iPlot));
        
        cfg = [];
        cfg.xlim = [time(iTimePoint),time(iTimePoint)];
        cfg.zlim = [minV,maxV];
        cfg.layout = 'elec1010.lay';
        cfg.comment = 'no';
        cfg.colormap = redblue;
        ft_topoplotER(cfg,data);
        
        % if iPlot <= nColSubplots, title(titleStr); end
        
        if mod(iPlot,nColSubplots) == 1 
            text(-0.1,0.5,strrep(conds{ceil(iPlot/nColSubplots)},'_','\_'),...
                'Units','normalized','HorizontalAlignment','right');
        end
        
        if iPlot == nSubplots
            switch orientation
                case 'horiz'
                    temp = findobj(h{iPlot},'Type','axes');
                    hAxes = temp(1);
                    posOutAxes = get(hAxes,'OuterPosition');
                    posOutCbar = [posOutAxes(1)+(posOutAxes(3)),posOutAxes(2),...
                        (posOutAxes(3)*0.2),posOutAxes(4)];
                    colorbar('Position',posOutCbar);
                case 'vert', colorbar('East');
            end
        end
        
        
        % switch orientation
            % case 'horiz', title(titleStr);
            % case 'vert', ylabel(titleStr);
        % end
    end
    
%     if ~isempty(weights)
%         for iPlot = 1:numel(h)-1, caxis(h{iPlot},[minV maxV]), end;
%         caxis(h{end},[minW maxW])
%     else
%         for iPlot = 1:numel(h), caxis(h{iPlot},[minV maxV]), end;
%     end
    for iPlot = 1:numel(h), caxis(h{iPlot},[minV maxV]), end;
    suplabel(sprintf('%d ms',round(time(iTimePoint)*1000)),'x',[.08 .12 .84 .84]);
    
end

end