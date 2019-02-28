function stats = compAVmodelCorrStats(indivData,time,genTime,varargin)
% Method to compute statistics for generalization performance.

% Parsing input
p = inputParser;

validGenTimes = {'tr','tr_x_tr'};

addRequired(p,'indivData',@(x) validateattributes(x,{'struct'},{'vector',...
    'nonempty'}));
addRequired(p,'time',@(x) validateattributes(x,{'numeric'},{'vector',...
    'increasing'}));
addRequired(p,'genTime',@(x) any(validatestring(x,validGenTimes)));
addParameter(p,'timeWin',[],@(x) validateattributes(x,{'numeric'},{'vector',...
    'increasing','numel',2}));
parse(p,indivData,time,genTime,varargin{:});

indivData = p.Results.indivData;
time = p.Results.time;
genTime = p.Results.genTime;
timeWin = p.Results.timeWin;

if isrow(time), time = time'; end

if isempty(timeWin)
    isInTimeWin = true(size(time));
else
    isInTimeWin = time >= timeWin(1) & time <= timeWin(2);
end

% Number of subjects
nSubj = length(indivData);
% Analysis names. As for now statistics are implemented only for main 
% effects
anNames = {'cc_acrossTime','cc_behav'};
anVarNames = {'cc_acrossTime','cc_behav'};

% Checking if any vairables are missing
dataFieldNames = fieldnames(indivData);
[isIn] = ismember(anVarNames,dataFieldNames);
anNames = anNames(isIn);
anVarNames = anVarNames(isIn);
% Array for collecting stats
tempStats = cell(size(anNames));
isVectorData = true(size(anNames));
% Full factorial expanson of the factor levels always follows this patern 
ff = [1,2];
nCond = size(ff,2);
indivData = squeeze(struct2cell(indivData));
for i = 1:numel(anNames)
    % Choosing variables for the actual analysis
    [~,actVarIdx] = ismember(anVarNames{i},dataFieldNames);
    temp = indivData(actVarIdx,:);
    % In case of a one sample test generate matching zero data
    if size(temp,1) == 1
        temp(2,:) = cellfun(@(x) zeros(size(x)),temp,'UniformOutput',false);
    end
    temp = temp(:)';
    isVectorData(i) = isvector(temp{1});
    data = cell(size(temp));
    % I'm using fieldtrip's high level statistical functions to perform the
    % cluster based permutation testing. Therefore the data must be 
    % converted to fieldtrip timelock or freq datasturcture. See the 
    % documentation of ft_datatype_timelock and ft_datatype_freq 
    % for details.
    if isVectorData(i)
        for j = 1:size(temp,2)
            data{j}.dimord = 'chan_time';
            data{j}.avg = temp{j}(isInTimeWin);
            if iscolumn(data{j}.avg), data{j}.avg = data{j}.avg'; end
            data{j}.label = {'foo'};
            data{j}.time = time(isInTimeWin);
        end
    else
        % In case data is generalized across time
        for j = 1:size(temp,2)
            data{j}.dimord = 'chan_freq_time';
            % Pretending there is 1 channel, T frequencies and T
            % timepoints, where T is the number of training time points
            data{j}.powspctrm = shiftdim(temp{j},-1);
            data{j}.label = {'foo'};
            data{j}.freq = 1:numel(time);
            data{j}.time = time;
        end
    end
    % Creating the design
    ivars = repmat(ff,1,nSubj);
    uvar = repmat(1:nSubj,nCond,1);
    uvar = uvar(:)';
    design = cat(1,uvar,ivars);
    % Configuration structure
    % settings for ft_timelockstatistics
    cfg.method           = 'montecarlo';
    cfg.avgoverchan      = 'yes';
    % settings for ft_statistics_montecarlo
    cfg.design           = design;
    cfg.numrandomization = 5000;
    cfg.correctm         = 'cluster';
    cfg.alpha            = 0.05;
    cfg.tail             = 0;
    cfg.ivar             = 2;
    cfg.uvar             = 1;
    cfg.randomseed       = 'yes';
    cfg.feedback         = 'text';
    cfg.clusterstatistic = 'maxsum';
    cfg.clusterthreshold = 'parametric';
    cfg.clusteralpha     = 0.05;
    cfg.clustertail      = 0;
    if cfg.tail == 0
        cfg.correcttail  = 'alpha';
    end
    cfg.statistic        = 'ft_statfun_depsamplesT';
    % Computing statistics 
    if isVectorData(i)
        cfg.parameter = 'avg';
        tempStats{i} = ft_timelockstatistics(cfg,data{:});
    else
        cfg.parameter = 'powspctrm';
        tempStats{i} = ft_freqstatistics(cfg,data{:});
    end
end

for i = 1:numel(anNames)
    pName = sprintf('p_%s_%s',genTime,anNames{i});
    hName = sprintf('h_%s_%s',genTime,anNames{i});
    stName = sprintf('st_%s_%s',genTime,anNames{i});
    if isVectorData(i)
        stats.(pName) = NaN(size(time));
        stats.(pName)(isInTimeWin) = tempStats{i}.prob';
        stats.(hName) = NaN(size(time));
        stats.(hName)(isInTimeWin) = tempStats{i}.mask';
        stats.(stName) = NaN(size(time));
        stats.(stName)(isInTimeWin) = tempStats{i}.stat';
    else
        stats.(pName) = squeeze(tempStats{i}.prob);
        stats.(hName) = squeeze(tempStats{i}.mask);
        stats.(stName) = squeeze(tempStats{i}.stat);
    end
end

end

% for i = 1:numel(anNames)
%     % Choosing variables for the actual analysis
%     [~,actVarIdx] = ismember(anVarNames{i},dataFieldNames);
%     temp = indivData(actVarIdx,:);
%     data = cat(2,temp{:});
%     [tempStats{i}.prob,tempStats{i}.stat,tempStats{i}.h] = deal(NaN(size(time)));
%     % Computing statistics
%     [~,~,pc,hc,st] = mvpa.bootstrpOneSampleTtest(data,10000,'MCPsol','cluster',...
%                                                  'clusterStat','maxsize');
%     tempStats{i}.prob = pc;
%     tempStats{i}.h = hc;
%     tempStats{i}.stat = st;

% end
% genTimeStr = 'tr';
% for i = 1:numel(anNames)
%     pName = sprintf('p_%s_%s',genTimeStr,anNames{i});
%     hName = sprintf('h_%s_%s',genTimeStr,anNames{i});
%     stName = sprintf('st_%s_%s',genTimeStr,anNames{i});
%     stats.(pName) = tempStats{i}.prob;
%     stats.(hName) = stats.(pName) < 0.05;
%     stats.(stName) = tempStats{i}.stat;
% end

% end


