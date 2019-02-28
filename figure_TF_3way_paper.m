function figure_TF_3way(varargin)
validPlotModes = {'channel','topo'};
p = inputParser;
addOptional(p,'plotMode','channel',@(x) ismember(x,validPlotModes));
parse(p,varargin{:})
plotMode = p.Results.plotMode;

expStage = 'final';
sourceDf = DEC_2_setupdir(expStage,'anal_eeg_group_tf');
saveDf = cd(sourceDf);

% Loading statistics
temp = load('stat_TF-pow-low_3way.mat');
stats = struct();
stats.low = temp.stats;
temp = load('stat_TF-pow-high_3way.mat');
stats.high = temp.stats;

% Loading TF data and plotting
conds = struct();
freqBands = {'low','high'};
for iFreq = 1:numel(freqBands)
    % Loading the 16 conditions of the 2x2x2x2 design
    listing = dir;
    listing = regexp({listing.name}',...
        ['fteeg_TF-pow-',freqBands{iFreq},'_',...
        'bis-r[hl]-d[hl]-[av]','_gravg.mat'],'match','once');
    fNames = listing(~strcmp(listing,''));
    listing = dir;
    tokens = regexp({listing.name}',...
        ['fteeg_TF-pow-',freqBands{iFreq},'_',...
        'bis-(r[hl]-d[hl]-[av])','_gravg.mat'],'tokens','once');
    tokens = tokens(~cellfun(@isempty,tokens));
    tokens = [tokens{:}]';
    temp = struct2cell(cellfun(@load,fNames));
        
    conds.(freqBands{iFreq}) = cat(2,strrep(tokens,'-','_'),temp');
    % Main effects
    matchStr = {...
        'VR',{'rh','rl'};...
        'Task',{'a','v'};...
        'Disp',{'dh','dl'}};
    c = cell(3,size(matchStr,1));
    for i = 1:size(matchStr,1)
        for j = 1:size(matchStr{i,2},2)
            idx = ~cellfun(@isempty,regexp(conds.(freqBands{iFreq})(:,1),...
                matchStr{i,2}{j},'match','once'));
            c{j,i} = avgFTdata(conds.(freqBands{iFreq}){idx,2});
        end
        c{3,i} = c{1,i};
        c{3,i}.powspctrm = c{2,i}.powspctrm-c{1,i}.powspctrm;
        % Averaging over subjects
        c(:,i) = cellfun(@(x) ft_freqdescriptives([],x),c(:,i),'UniformOutput',false);
        % Baseline correcting main effect levels
        cfg = struct();
        cfg.baseline = [-0.4,-0.2];
        c(1:2,i) = cellfun(@(x) ft_freqbaseline(cfg,x),c(1:2,i),'UniformOutput',false);
        % Adding mask from stats
        c{1,i}.mask = false(size(stats.(freqBands{iFreq}).(matchStr{i,1}).mask));
        c{2,i}.mask = false(size(stats.(freqBands{iFreq}).(matchStr{i,1}).mask));
        c{3,i}.mask = stats.(freqBands{iFreq}).(matchStr{i,1}).mask;
    end
    main.(freqBands{iFreq}) = array2table(c,'VariableNames',matchStr(:,1),...
        'RowNames',{'A1','A2','eff'});
    
    % 2-way interactions
    matchStr = {...
        'VR_X_Task',{'rh_.*_a','rh_.*_v','rl_.*_a','rl_.*_v'};...
        'VR_X_Disp',{'rh_dh','rh_dl','rl_dh','rl_dl'};...
        'Disp_X_Task',{'dh_a','dh_v','dl_a','dl_v'}};
    c = cell(9,size(matchStr,1));
    for i = 1:size(matchStr,1)
        for j = 1:size(matchStr{i,2},2)
            idx = ~cellfun(@isempty,regexp(conds.(freqBands{iFreq})(:,1),...
                matchStr{i,2}{j},'match','once'));
            c{j,i} = avgFTdata(conds.(freqBands{iFreq}){idx,2});
        end
        c{5,i} = c{1,i};
        c{5,i}.powspctrm = c{3,i}.powspctrm-c{1,i}.powspctrm;
        c{6,i} = c{2,i};
        c{6,i}.powspctrm = c{4,i}.powspctrm-c{2,i}.powspctrm;
        c{7,i} = c{1,i};
        c{7,i}.powspctrm = c{2,i}.powspctrm-c{1,i}.powspctrm;
        c{8,i} = c{3,i};
        c{8,i}.powspctrm = c{4,i}.powspctrm-c{3,i}.powspctrm;
        c{9,i} = c{1,i};
        c{9,i}.powspctrm = (c{4,i}.powspctrm-c{3,i}.powspctrm) - ...
            (c{2,i}.powspctrm-c{1,i}.powspctrm);
        % Averaging over subjects
        c(:,i) = cellfun(@(x) ft_freqdescriptives([],x),c(:,i),'UniformOutput',false);
        % Adding mask from stats
        for k = 1:size(c,1)-1
            c{k,i}.mask = false(size(stats.(freqBands{iFreq}).(matchStr{i,1}).mask));
        end
        c{9,i}.mask = stats.(freqBands{iFreq}).(matchStr{i,1}).mask;
    end
    int2way.(freqBands{iFreq}) = array2table(c,'VariableNames',matchStr(:,1),...
        'RowNames',{'A1B1','A1B2','A2B1','A2B2','AinB1','AinB2','BinA1',...
                    'BinA2','eff'});
    
end
cd(saveDf);

% % Plotting figures
if strcmp(plotMode,'channel')
    chanClusters = {...
        {'FC5','FC3','FC1','FC2','FC4','FC6','F5','F3','F1','Fz','F2','F4','F6'},...
        {'P5','P3','P1','Pz','P2','P4','P6','CP5','CP3','CPz','CP1','CP2',...
        'CP4','CP6','C5','C3','C1','Cz','C2','C4','C6'},...
        {'O1','Oz','O2','PO7','PO3','POz','PO4','PO8'}};
    chanClusters = {...
        {'F1','F3','F5','FC1','FC3','FC5','F2','F4','F6','FC2','FC4','FC6'};...
        {'C1','C3','C5','CP1','CP3','CP5','C2','C4','C6','CP2','CP4','CP6'};...
        {'O1','PO3','PO7','P1','P3','P5','O2','PO4','PO8','P2','P4','P6'}}';
    chanClusters = {{'Fz'},{'Cz'},{'Pz'},{'Oz'}};
    clusterNames = {'Frontal','Central','Parietal','Occipital'};
    % freqRanges = {[4,8],[8,12],[13,30]};
    freqRanges = {[8,30]};
    effs = {'VR','Task','Disp','VR_X_Task'};
    figure();
    for iEff = 1:numel(effs)
        for i = 1:numel(freqRanges)
            ax = {};
            for j = 1:numel(chanClusters)
                cfg = struct();
                cfg.frequency = freqRanges{i};
                cfg.channel = chanClusters{j}';
                cfg.avgoverfreq = 'yes';
                cfg.nanmean = 'yes';
                temp = {};
                if iEff < 4
                    for k = 1:3
                        temp{k} = ft_selectdata(cfg,main.low.(effs{iEff}){k});
                    end
                else
                    rowsToPlot = {'AinB1','AinB2','eff'};
                    for k = 1:3
                        temp{k} = ft_selectdata(cfg,int2way.low.(effs{iEff}){rowsToPlot{k}});
                    end
                end
%                 ax{j} = subplot(numel(effs),6,(6*(iEff-1))+(3*(i-1))+j);
                ax{j} = subplot(numel(effs),numel(chanClusters),...
                    (numel(chanClusters)*(iEff-1))+j);
                cfg = struct();
    
                cfg.ylim = 'maxabs';
                cfg.xlim = [-0.4,0.7];
                cfg.baseline = 'no';
                cfg.maskparameter = 'mask';
                cfg.graphcolor = 'grb';
                ft_singleplotER(cfg,temp{numel(temp):-1:1});
                if i == 1 && j == 1
                    ylabel(strrep(effs{iEff},'_',' '));
                end
                if iEff > 1
                    title('');
                end
                if j > 1
                    set(gca,'YTick',[]);
                end
            end
            ax = cat(2,ax{:});
            for iAx = 1:size(ax,2)
                setYlimits(ax);
            end
        end
    end
    
else
    % Clusterplots
    % Flip the sign of statistics because the contrast is computed the other
    % way around here
    effs = {'VR','Task','Disp'};
    for i = 1:numel(effs)
        temp = stats.low.(effs{i}).stat;
        stats.low.(effs{i}).stat = -temp;
    end
    
    effs = {'VR','Task','Task','VR_X_Task'};
    freqRanges = {[8,30],[8,30],[8,30],[8,30]};
    freqStr = {'alpha-beta','alpha-beta','alpha-beta','alpha-beta'};
    latencies = {[0.2,0.4],[-0.2,0],[0.35,0.55],[0.4,0.6]};
    figure();
    colormap('redblue');
    % colormap('parula');
    plotOrder = [1,3,4,5];
    ax = {};
    for i = 1:numel(effs)
        cfg = struct();
        cfg.frequency = freqRanges{i};
        cfg.latency = 'all';
        cfg.avgoverfreq = 'yes';
        cfg.latency = latencies{i};
        cfg.avgovertime = 'yes';
        cfg.nanmean = 'yes';
        temp = ft_selectdata(cfg,stats.low.(effs{i}));
        % Selecting statistics within the time window of interest. If a
        % channel has a significant value somewhere within the time window,
        % the channel will be significant.
        cfg = struct();
        cfg.frequency = freqRanges{i};
        cfg.latency = latencies{i};
        tempStat = ft_selectdata(cfg,stats.low.(effs{i}));
        tempMask = any(any(tempStat.mask,3),2);
        temp.mask = tempMask;
        
        cfg = struct();
        cfg.layout = 'acticap-64ch-standard2';
        cfg.parameter = 'stat';
        cfg.colorbar = 'no';
        cfg.zlim = 'maxabs';
        cfg.comment = 'no';
        cfg.highlight = 'on';
        cfg.highlightchannel = temp.label(temp.mask);
        cfg.highlightsymbol = '*';
        cfg.highlightsize = 4;
        ax{i} = subplot(3,2,plotOrder(i));
        ft_topoplotTFR(cfg,temp);
        title(sprintf('%s band\n%d - %d ms',...
            freqStr{i},...
            latencies{i}(1)*1000,latencies{i}(2)*1000));
        if mod(plotOrder(i),2) == 1
            text(-1,0,strrep(effs{i},'_',' '),'Rotation',90,...
                'HorizontalAlignment','center');
        end
    end
    c = colorbar('EastOutside');
    c.Label.String = 't-value';
    % Setting color
    ax = cat(2,ax{:});
    for iAx = 1:size(ax,2)
        setColorLimits(ax);
    end
end

end

function ylim = setYlimits(ax,varargin)
if isempty(varargin)
    temp = arrayfun(@(i) get(ax(i),'ylim'),1:numel(ax),'UniformOutput',false);
    temp = cat(2,temp{:});
    ylim = max(abs(min(temp)),abs(max(temp)));
else
    ylim = varargin{1};
end
arrayfun(@(i) set(ax(i),'ylim',[-ylim,ylim]),1:numel(ax));
end

function clim = setColorLimits(ax,varargin)
if isempty(varargin)
    temp = arrayfun(@(i) get(ax(i),'clim'),1:numel(ax),'UniformOutput',false);
    temp = cat(2,temp{:});
    clim = max(abs(min(temp)),abs(max(temp)));
else
    clim = varargin{1};
end
arrayfun(@(i) set(ax(i),'clim',[-clim,clim]),1:numel(ax));
end