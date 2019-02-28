function stats = compEEGstats(indivData,expDesign,dataType)
% Function to compute statistics for time-frequency (TF) and ERP data.

% Parsing input
validExpDesigns = {'oneSample','3way','4way'};
validDataTypes = {'freq','timelock'};
p = inputParser;
addRequired(p,'indivData',@(x) validateattributes(x,{'struct'},{'vector',...
    'nonempty'}));
addRequired(p,'expDesign',@(x) ismember(x,validExpDesigns));
addRequired(p,'dataType',@(x) ismember(x,validDataTypes));
parse(p,indivData,expDesign,dataType);
indivData = p.Results.indivData;
expDesign = p.Results.expDesign;
dataType = p.Results.dataType;

% Number of subjects
dataFieldNames = fieldnames(indivData);
% Assuming indiv data was created from ft_freqgrandaverage and each
% condition has the same number of subjects
if strcmp(dataType,'freq')
    nSubj = size(indivData.(dataFieldNames{1}).powspctrm,1);
else
    nSubj = size(indivData.(dataFieldNames{1}).individual,1);
end
if strcmp(expDesign,'oneSample')
    anEffects = {'Overall'};                      % 3-way ineraction
    % I perform statistics on powerspectra from the 8 conditions of the
    % 3-way design
    anVars = {{'rh_dh_a','rh_dh_v','rh_dl_a','rh_dl_v',...
              'rl_dh_a','rl_dh_v','rl_dl_a','rl_dl_v'}};
    % Standardized names for effects - not necessary here
    anEffectsStd = {'bla'};
    % Variable names converted to standardized factor names - not necessary
    % here
    anVarsStd = {{'bla'}};
elseif strcmp(expDesign,'3way')
    anEffects = {'VR','Task','Disp',...                    % main effects
                 'VR_X_Disp','VR_X_Task','Disp_X_Task',... % 2-way interactions
                 'VR_X_Task_X_Disp'};                      % 3-way ineraction
    % I perform statistics on powerspectra from the 8 conditions of the
    % 3-way design
    anVars = repmat({{'rh_dh_a','rh_dh_v','rh_dl_a','rh_dl_v',...
                      'rl_dh_a','rl_dh_v','rl_dl_a','rl_dl_v'}},1,7);
    % Standardized names for effects
    anEffectsStd = {'A','C','B','A X B','A X C','B X C','A X B X C'}; % 3-way design
    
    % Variable names converted to standardized factor names
    anVarsStd = repmat({{'A1_B1_C1','A1_B1_C2','A1_B2_C1','A1_B2_C2',...
                         'A2_B1_C1','A2_B1_C2','A2_B2_C1','A2_B2_C2'}},1,7);
else
    anEffects = {'VR','Task','VLvsR','ALvsR',...    % main effects 
        'VR_X_Task','VR_X_VLvsR','VR_X_ALvsR',...   % 2-way int 
        'Task_X_VLvsR','Task_X_ALvsR'};
    % I perform statistics on powerspectra from the 16 conditions of 
    % the 4-way design
    anVars = repmat(...
        {{'rh_a_vl_al','rh_a_vl_ar','rh_a_vr_al','rh_a_vr_ar',...
        'rh_v_vl_al','rh_v_vl_ar','rh_v_vr_al','rh_v_vr_ar',...
        'rl_a_vl_al','rl_a_vl_ar','rl_a_vr_al','rl_a_vr_ar',...
        'rl_v_vl_al','rl_v_vl_ar','rl_v_vr_al','rl_v_vr_ar'}},1,9);
    % Standardized names for effects
    anEffectsStd = {'A','B','C','D','A X B','A X C','A X D','B X C','B X D'};
    % Variable names converted to standardized factor names
    anVarsStd = repmat(...
        {{'A1_B1_C1_D1','A1_B1_C1_D2','A1_B1_C2_D1','A1_B1_C2_D2',...
        'A1_B2_C1_D1','A1_B2_C1_D2','A1_B2_C2_D1','A1_B2_C2_D2',...
        'A2_B1_C1_D1','A2_B1_C1_D2','A2_B1_C2_D1','A2_B1_C2_D2',...
        'A2_B2_C1_D1','A2_B2_C1_D2','A2_B2_C2_D1','A2_B2_C2_D2'}},1,9);
end

tempStats = cell(size(anEffects));
indivData = struct2cell(indivData);
% prepare neighbouring channel groups
cfg = struct();
cfg.feedback = 'no';
cfg.method = 'triangulation';
cfg.layout = 'acticap-64ch-standard2';
neighbours = ft_prepare_neighbours(cfg,indivData{1});
parfor i = 1:numel(anEffects)
    % Choosing variables for the actual analysis
    [~,actVarIdx] = ismember(anVars{i},dataFieldNames);
    data = indivData(actVarIdx,:);
    % Configuration structure
    cfg = struct();
    % settings for ft_freqstatistics
    cfg.method = 'montecarlo'; 
    % settings for ft_statistics_montecarlo
    % Creating the design
    if strcmp(expDesign,'oneSample')
        [cfg,data] = buildDesignOneSample(cfg,data,nSubj);
    elseif strcmp(expDesign,'3way')
        [cfg,data] = buildDesignThreeWay(cfg,data,nSubj,anEffectsStd{i},anVarsStd{i});
    elseif strcmp(expDesign,'4way')
        [cfg,data] = buildDesignFourWay(cfg,data,nSubj,anEffectsStd{i},anVarsStd{i});
    end
    cfg.numrandomization = 'all';
    cfg.correctm = 'cluster'; %'max'
    cfg.clusterstatistic = 'maxsum';
    cfg.clusterthreshold = 'parametric';
    cfg.clusteralpha = 0.05;
    cfg.clustertail = 0;
    cfg.alpha = 0.05;
    cfg.tail = 0;
    cfg.correcttail = 'alpha';
    cfg.randomseed = 'yes';
    cfg.feedback = 'text';
    cfg.statistic = 'ft_statfun_depsamplesT';
    cfg.neighbours = neighbours;
    % Computing statistics 
    if strcmp(dataType,'freq')
        cfg.parameter = 'powspctrm';
        tempStats{i} = ft_freqstatistics(cfg,data{:});
    else
        cfg.parameter = 'individual';
        tempStats{i} = ft_timelockstatistics(cfg,data{:});
    end
end

stats = cell2struct(tempStats,anEffects,2);

end

% Accessory function for building the one sample design
function [cfg,data] = buildDesignOneSample(cfg,data,nSubj)
 
% Averaging over all conditions
data_avg = avgFTdata(data{:});

if isfield(data_avg,'individual')
    parameter = 'individual';
else
    parameter = data_avg.cfg.parameter{:};
end
data_zero = data_avg;
temp = data_zero.(parameter);
temp(~isnan(temp)) = 0;
data_zero.(parameter) = temp;
data = {data_avg,data_zero};
% Independent variable
ivar = repmat([1,2],nSubj,1);
ivar = ivar(:)';
% Unit variable
uvar = repmat((1:nSubj)',1,numel(data));
uvar = uvar(:)';
% Version with averaging
cfg.design = cat(1,uvar,ivar);
cfg.uvar = 1;
cfg.ivar = 2;

end

% Accessory function for building the 3-way design
function [cfg,data] = buildDesignThreeWay(cfg,data,nSubj,anName,varNames)
% varNames must be defined in order A_B_C
nCond = size(varNames,2);
% Split variable names to get individual factor levels
varConds = cellfun(@strsplit,varNames,repmat({'_'},size(varNames)),'UniformOutput',false);
% nConditions x nFactors
varConds = cat(1,varConds{:});
 
if ismember(anName,{'A','B','C'})
%     % Independent variables. One row for each factor, recoded to [1,2]
%     % within each level
%     ivar = ones(size(varConds));
%     ivar(ismember(varConds,{'A2','B2','C2'})) = 2;
%     ivar = repmat(ivar,1,1,nSubj);
%     ivar = shiftdim(ivar,1);
%     ivar = reshape(ivar,3,[]);
%     % Unit variable. This indicates the unit of observation, the
%     % subject. Values are 1:nSubj. This restricts permutations within
%     % each subject.
%     uvar = repmat((1:nSubj)',1,nCond);
%     uvar = uvar(:)';
    % The second dimension restricts the permutations to within certain
    % conditions based on the analysis at hand.
    % Main effects. Permute separetely within the levels of the not
    % tested effects. 
    if strcmp(anName,'A')
        % Test main effect of A, permute within levels of B and C.
        idx = ismember(varConds(:,1),'A1')';
        varIDx = 1;
    elseif strcmp(anName,'B')
        % Test main effect of B, permute within levels of A and C.
        idx = ismember(varConds(:,2),'B1')';
        varIDx = 2;
    elseif strcmp(anName,'C')
        % Test main effect of C, permute within levels of A and B.
        idx = ismember(varConds(:,3),'C1')';
        varIDx = 3;
    end
    % Averaging over not-tested conditions
    data_avgs{1} = avgFTdata(data{idx});
    data_avgs{2} = avgFTdata(data{~idx});
    data = data_avgs;
    % Independent variable
    ivar = repmat([1,2],nSubj,1);
    ivar = ivar(:)';
    % Unit variable
    uvar = repmat((1:nSubj)',1,numel(data));
    uvar = uvar(:)';
%     % Within-cell variable: examples within the levels of this variable are
%     % kept together during permutations
%     wvar = ivar(varIDx,:);
    % Saving design and related info. I kick out other factors from ivar
    % just for good measure, as depsamplesT expects one factor. 
    % Version with wvar
%     cfg.design = cat(1,uvar,wvar,ivar(varIDx,:));
%     cfg.uvar = 1;
%     cfg.wvar = 2;
%     cfg.ivar = 3;
% Version with averaging
    cfg.design = cat(1,uvar,ivar);
    cfg.uvar = 1;
    cfg.ivar = 2;
elseif ismember(anName,{'A X B','A X C','B X C'}) 
    % Two-way interactions. Permute the differences of the levels of the
    % two factors (11-12),(21-22) within levels of 3rd factors
    % Conditions organized into differences, which are taken from each row
    if strcmp(anName,'A X B')
        diffs = {'A1_B1_C1','A1_B2_C1';...
                 'A1_B1_C2','A1_B2_C2';...
                 'A2_B1_C1','A2_B2_C1';...
                 'A2_B1_C2','A2_B2_C2'};
    elseif strcmp(anName,'A X C')
        diffs = {'A1_B1_C1','A1_B1_C2';...
                 'A1_B2_C1','A1_B2_C2';...
                 'A2_B1_C1','A2_B1_C2';...
                 'A2_B2_C1','A2_B2_C2'};
    elseif strcmp(anName,'B X C')
        diffs = {'A1_B1_C1','A1_B1_C2';...
                 'A2_B1_C1','A2_B1_C2';...
                 'A1_B2_C1','A1_B2_C2';...
                 'A2_B2_C1','A2_B2_C2'};
    end
    % Creating difference datasets following fieldtrip FAQ
    % https://tinyurl.com/yagc5yey
    data_diffs = cell(size(diffs,1),1);
    for i = 1:size(diffs,1)
        data_diffs{i} = diffFTdata(data{ismember(varNames,diffs{i,1})},...
                                   data{ismember(varNames,diffs{i,2})});
    end
%     data = data_diffs;
    % Averaging over not-tested conditions
    data_avgs{1} = avgFTdata(data_diffs{1:2});
    data_avgs{2} = avgFTdata(data_diffs{3:4});
    data = data_avgs;
    % Permutation from this point is the same as for main effects: permute
    % the difference datasets within the factors of the third level
    % Independent variables. 
%     ivar = repmat([1,1,2,2],nSubj,1);
    ivar = repmat([1,2],nSubj,1);
    ivar = ivar(:)';
    % Unit variable. 
    uvar = repmat((1:nSubj)',1,numel(data));
    uvar = uvar(:)';
%     % Within-cell variable
%     wvar = ivar;
    % Saving design and related info. 
    % Version with wvar
%     cfg.design = cat(1,uvar,wvar,ivar);
%     cfg.uvar = 1;
%     cfg.wvar = 2;
%     cfg.ivar = 3;
    % Version with averaging
    cfg.design = cat(1,uvar,ivar);
    cfg.uvar = 1;
    cfg.ivar = 2;
elseif strcmp(anName,'A X B X C')
    % Three-way interaction, permute freely
    % Computing one two-way interaction
    diffs = {'A1_B1_C1','A1_B2_C1';...
             'A1_B1_C2','A1_B2_C2';...
             'A2_B1_C1','A2_B2_C1';...
             'A2_B1_C2','A2_B2_C2'};
    data_diffs = cell(size(diffs,1),1);
    for i = 1:size(diffs,1)
        data_diffs{i} = diffFTdata(data{ismember(varNames,diffs{i,1})},...
                                    data{ismember(varNames,diffs{i,2})});
    end
    % Computing the three-way interaction from that
    temp{1} = diffFTdata(data_diffs{1},data_diffs{3});
    temp{2} = diffFTdata(data_diffs{2},data_diffs{4});
    data = temp;
    % Independent variables. 
    ivar = repmat([1,2],nSubj,1);
    ivar = ivar(:)';
    % Unit variable. 
    uvar = repmat((1:nSubj)',1,numel(data));
    uvar = uvar(:)';
    % Saving design and related info
    cfg.design = cat(1,uvar,ivar);
    cfg.uvar = 1;
    cfg.ivar = 2;
    
else
    error('compTFstats:unknownEffect',...
          'This effect is not implemented!')
end

end

% Accessory function for building the 4-way design
function [cfg,data] = buildDesignFourWay(cfg,data,nSubj,anName,varNames)
% varNames must be defined in order A_B_C_D
nCond = size(varNames,2);
% Split variable names to get individual factor levels
varConds = cellfun(@strsplit,varNames,repmat({'_'},size(varNames)),'UniformOutput',false);
% nConditions x nFactors
varConds = cat(1,varConds{:});
 
if ismember(anName,{'A','B','C','D'})
%     % Independent variables. One row for each factor, recoded to [1,2]
%     % within each level
%     ivar = ones(size(varConds));
%     ivar(ismember(varConds,{'A2','B2','C2','D2'})) = 2;
%     ivar = repmat(ivar,1,1,nSubj);
%     ivar = shiftdim(ivar,1);
%     ivar = reshape(ivar,4,[]);
%     % Unit variable. This indicates the unit of observation, the
%     % subject. Values are 1:nSubj. This restricts permutations within
%     % each subject.
%     uvar = repmat((1:nSubj)',1,nCond);
%     uvar = uvar(:)';
    % Main effects. Permute separetely within the levels of the not
    % tested effects. 
    if strcmp(anName,'A')
        % Test main effect of A, permute within levels of B, C and D.
        idx = ismember(varConds(:,1),'A1')';
        varIDx = 1;
    elseif strcmp(anName,'B')
        % Test main effect of B, permute within levels of A, C and D.
        idx = ismember(varConds(:,2),'B1')';
        varIDx = 2;
    elseif strcmp(anName,'C')
        % Test main effect of C, permute within levels of A, B and D.
        idx = ismember(varConds(:,3),'C1')';
        varIDx = 3;
    elseif strcmp(anName,'D')
        % Test main effect of C, permute within levels of A, B and C.
        idx = ismember(varConds(:,4),'D1')';
        varIDx = 4;
    end
    % Averaging over not-tested conditions
    data_avgs{1} = avgFTdata(data{idx});
    data_avgs{2} = avgFTdata(data{~idx});
    data = data_avgs;
    % Independent variable
    ivar = repmat([1,2],nSubj,1);
    ivar = ivar(:)';
    % Unit variable
    uvar = repmat((1:nSubj)',1,numel(data));
    uvar = uvar(:)';
%     % Within-cell variable: examples within the levels of this variable are
%     % kept together during permutations
%     wvar = ivar(varIDx,:);
    % Saving design and related info. I kick out other factors from ivar
    % just for good measure, as depsamplesT expects one factor. 
% Version with wvar
%     cfg.design = cat(1,uvar,wvar,ivar(varIDx,:));
%     cfg.uvar = 1;
%     cfg.wvar = 2;
%     cfg.ivar = 3;
% Version with averaging
    cfg.design = cat(1,uvar,ivar);
    cfg.uvar = 1;
    cfg.ivar = 2;
elseif ismember(anName,{'A X B','A X C','A X D','B X C','B X D'}) 
    % Two-way interactions. Permute the differences of the levels of the
    % two factors (11-12),(21-22) within levels of 3rd factors
    % Conditions organized into differences, which are taken from each row
    if strcmp(anName,'A X B')
        diffs = {'A1_B1_C1_D1','A1_B2_C1_D1';...
                 'A1_B1_C1_D2','A1_B2_C1_D2';...
                 'A1_B1_C2_D1','A1_B2_C2_D1';...
                 'A1_B1_C2_D2','A1_B2_C2_D2';...
                 'A2_B1_C1_D1','A2_B2_C1_D1';...
                 'A2_B1_C1_D2','A2_B2_C1_D2';...
                 'A2_B1_C2_D1','A2_B2_C2_D1';...
                 'A2_B1_C2_D2','A2_B2_C2_D2'};
    elseif strcmp(anName,'A X C')
        diffs = {'A1_B1_C1_D1','A1_B1_C2_D1';...
                 'A1_B1_C1_D2','A1_B1_C2_D2';...
                 'A1_B2_C1_D1','A1_B2_C2_D1';...
                 'A1_B2_C1_D2','A1_B2_C2_D2';...
                 'A2_B1_C1_D1','A2_B1_C2_D1';...
                 'A2_B1_C1_D2','A2_B1_C2_D2';...
                 'A2_B2_C1_D1','A2_B2_C2_D1';...
                 'A2_B2_C1_D2','A2_B2_C2_D2'};
    elseif strcmp(anName,'A X D')
        diffs = {'A1_B1_C1_D1','A1_B1_C1_D2';...
                 'A1_B1_C2_D1','A1_B1_C2_D2';...
                 'A1_B2_C1_D1','A1_B2_C1_D2';...
                 'A1_B2_C2_D1','A1_B2_C2_D2';...
                 'A2_B1_C1_D1','A2_B1_C1_D2';...
                 'A2_B1_C2_D1','A2_B1_C2_D2';...
                 'A2_B2_C1_D1','A2_B2_C1_D2';...
                 'A2_B2_C2_D1','A2_B2_C2_D2'};
    elseif strcmp(anName,'B X C')
        diffs = {'A1_B1_C1_D1','A1_B1_C2_D1';...
                 'A1_B1_C1_D2','A1_B1_C2_D2';...
                 'A2_B1_C1_D1','A2_B1_C2_D1';...
                 'A2_B1_C1_D2','A2_B1_C2_D2';...
                 'A1_B2_C1_D1','A1_B2_C2_D1';...
                 'A1_B2_C1_D2','A1_B2_C2_D2';...
                 'A2_B2_C1_D1','A2_B2_C2_D1';...
                 'A2_B2_C1_D2','A2_B2_C2_D2'};
    elseif strcmp(anName,'B X D')
        diffs = {'A1_B1_C1_D1','A1_B1_C1_D2';...
                 'A1_B1_C2_D1','A1_B1_C2_D2';...
                 'A2_B1_C1_D1','A2_B1_C1_D2';...
                 'A2_B1_C2_D1','A2_B1_C2_D2';...
                 'A1_B2_C1_D1','A1_B2_C1_D2';...
                 'A1_B2_C2_D1','A1_B2_C2_D2';...
                 'A2_B2_C1_D1','A2_B2_C1_D2';...
                 'A2_B2_C2_D1','A2_B2_C2_D2'};
    end
    % Creating difference datasets following fieldtrip FAQ
    % https://tinyurl.com/yagc5yey
    data_diffs = cell(size(diffs,1),1);
    for i = 1:size(diffs,1)
        data_diffs{i} = diffFTdata(data{ismember(varNames,diffs{i,1})},...
                                   data{ismember(varNames,diffs{i,2})});
    end
%     data = data_diffs;
    % Averaging over not-tested conditions
    data_avgs{1} = avgFTdata(data_diffs{1:4});
    data_avgs{2} = avgFTdata(data_diffs{5:8});
    data = data_avgs;
    % Permutation from this point is the same as for main effects: permute
    % the difference datasets within the factors of the third level
    % Independent variables. 
%     ivar = repmat([1,1,1,1,2,2,2,2],nSubj,1);
    ivar = repmat([1,2],nSubj,1);
    ivar = ivar(:)';
    % Unit variable. 
    uvar = repmat((1:nSubj)',1,numel(data));
    uvar = uvar(:)';
%     % Within-cell variable
%     wvar = ivar;
    % Saving design and related info. 
%     % Version with wvar
%     cfg.design = cat(1,uvar,wvar,ivar);
%     cfg.uvar = 1;
%     cfg.wvar = 2;
%     cfg.ivar = 3;
    % Version with averaging
    cfg.design = cat(1,uvar,ivar);
    cfg.uvar = 1;
    cfg.ivar = 2;
else
    error('compTFstats:unknownEffect',...
          'This effect is not implemented!')
end

end