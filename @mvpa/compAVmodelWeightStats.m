function stats = compAVmodelWeightStats(indivData,time,varargin)
% Method to compute statistics for AV weights data.

% Parsing input
p = inputParser;

validModes = {'behav','eeg','eeg_sim'};
validAnalyses = {'VR X Task X Disp','Psalpha'};

addRequired(p,'indivData',@(x) validateattributes(x,{'struct'},{'vector',...
    'nonempty'}));
addRequired(p,'time',@(x) validateattributes(x,{'numeric'},{'vector',...
    'increasing'}));
addParameter(p,'mode','eeg',@(x) any(validatestring(x,validModes)));
addParameter(p,'analysis','VR X Task X Disp',...
             @(x) any(validatestring(x,validAnalyses)));
addParameter(p,'timeWin',[],@(x) validateattributes(x,{'numeric'},{'vector',...
    'increasing','numel',2}));
parse(p,indivData,time,varargin{:});

indivData = p.Results.indivData;
time = p.Results.time;
mode = p.Results.mode;
analysis = p.Results.analysis;
timeWin = p.Results.timeWin;

temp = fieldnames(indivData);
if length(time) ~= length(indivData(1).(temp{1}))
    error('mvpa:compAVmodelWeightStats:missingInput',...
          'Time must be the same length as the individual data');
end

if iscolumn(time), time = time'; end

if isempty(timeWin)
    isInTimeWin = true(size(time));
else
    isInTimeWin = time >= timeWin(1) & time <= timeWin(2);
end

% Number of subjects
nSubj = length(indivData);
% Analysis names.
if ismember(mode,{'behav','eeg'})
    if strcmp(analysis,'VR X Task X Disp')
        anEffects = {'VR','Task','Disp',...                 % main effects
                     'VR_X_Disp','VR_X_Task','Disp_X_Task',... % 2-way interactions
                     'VR_X_Task_X_Disp',...
                     'Task_in_DispLow','Task_in_DispHigh',...
                     'VR_in_A','VR_in_V',...
                     'VR_in_DispLow','VR_in_DispHigh'};
        % anEffects = {'Task_|_DispLow','Task_|_DispHigh',...
                     % 'VR_|_A','VR_|_V',...
                     % 'VR_|_DispLow','VR_|_DispHigh'};
        % For main effects and interactions we use all eight variables
        anVars = repmat({{'r_d_a_wav','r_d_v_wav','r_D_a_wav','r_D_v_wav',...
                         'R_d_a_wav','R_d_v_wav','R_D_a_wav','R_D_v_wav'}},...
                         1,7);
        % For simple effects we reduce the design into a 2 x 2
        % within the levels of the third factor
        % anVars = {{'r_d_a_wav','R_d_a_wav','r_d_v_wav','R_d_v_wav'}};
        anVars = cat(2,anVars,{{'r_d_a_wav','R_d_a_wav','r_d_v_wav','R_d_v_wav'}});
        anVars = cat(2,anVars,{{'r_D_a_wav','R_D_a_wav','r_D_v_wav','R_D_v_wav'}});
        anVars = cat(2,anVars,{{'r_d_a_wav','r_D_a_wav','R_d_a_wav','R_D_a_wav'}});
        anVars = cat(2,anVars,{{'r_d_v_wav','r_D_v_wav','R_d_v_wav','R_D_v_wav'}});
        anVars = cat(2,anVars,{{'r_d_a_wav','r_d_v_wav','R_d_a_wav','R_d_v_wav'}});
        anVars = cat(2,anVars,{{'r_D_a_wav','r_D_v_wav','R_D_a_wav','R_D_v_wav'}});
        
        % Standardized names for effects
        anEffectsStd = {'A','C','B','A X B','A X C','B X C','A X B X C',...
                        'A','A','A','A','A','A'};
        % anEffectsStd = {'A','A','A','A','A','A'};
        % Variable names converted to standardized factor names
        anVarsStd = repmat({{'A1_B1_C1','A1_B1_C2','A1_B2_C1','A1_B2_C2',...
                            'A2_B1_C1','A2_B1_C2','A2_B2_C1','A2_B2_C2'}},...
                           1,7);
        anVarsStd = cat(2,anVarsStd,...
                        repmat({{'A1_B1','A1_B2','A2_B1','A2_B2'}},1,6));
        % anVarsStd = repmat({{'A1_B1','A1_B2','A2_B1','A2_B2'}},1,6);
        expDesign = cat(2,repmat({'3way'},1,7),repmat({'2way'},1,6));
        % expDesign = repmat({'2way'},1,6);
    else
        anEffects = {'Psalpha'};
        anVars = {'psalow_wav','psahigh_wav'};
        % Standardized names for effects
        anEffectsStd = {'A'};
        % Variable names converted to standardized factor names
        anVarsStd = {'A1','A2'};
        expDesign = {'1way'};
    end
    
elseif strcmp(mode,'eeg_sim')
    % In simulation mode test only one for each of the three effect
    % types: main, 2-way and 3-way interactions
    anEffects = {'VR','VR_X_Task','VR_X_Task_X_Disp'};
    anEffectsStd = {'A','A X C','A X B X C'};
    
    anVars = repmat({{'r_d_a_wav','r_d_v_wav','r_D_a_wav','r_D_v_wav',...
                      'R_d_a_wav','R_d_v_wav','R_D_a_wav','R_D_v_wav'}},1,3);
    
    % Variable names converted to standardized factor names
    anVarsStd = repmat({{'A1_B1_C1','A1_B1_C2','A1_B2_C1','A1_B2_C2',...
                        'A2_B1_C1','A2_B1_C2','A2_B2_C1','A2_B2_C2'}},1,3);
    expDesign = repmat({'3way'},1,3);
end

tempStats = cell(size(anEffects));

dataFieldNames = fieldnames(indivData);
indivData = squeeze(struct2cell(indivData));
parfor i = 1:numel(anEffects)
    % Choosing variables for the actual analysis
    [~,actVarIdx] = ismember(anVars{i},dataFieldNames);
    temp = indivData(actVarIdx,:);
    temp = temp(:)';
    data = cell(size(temp));
    % I'm using fieldtrip's high level statistical functions to perform the
    % cluster based permutation testing. Therefore the data must be 
    % converted to fieldtrip timelock datasturcture. See the 
    % documentation of ft_datatype_timelock for details.
    for j = 1:size(temp,2)
        data{j}.dimord = 'chan_time';
        data{j}.avg = temp{j}(isInTimeWin);
        if iscolumn(data{j}.avg), data{j}.avg = data{j}.avg'; end
        data{j}.label = {'foo'};
        data{j}.time = time(isInTimeWin);
    end
    % Configuration structure
    cfg = struct();
    % settings for ft_timelockstatistics
    cfg.method = 'montecarlo'; 
    cfg.avgoverchan = 'yes';
    cfg.parameter = 'avg';
    % settings for ft_statistics_montecarlo
    % Creating the design
    if strcmp(expDesign{i},'3way')
        cfg = buildDesignThreeWay(cfg,nSubj,anEffectsStd{i},anVarsStd{i});
    elseif strcmp(expDesign{i},'2way')
        cfg = buildDesignTwoWay(cfg,nSubj,anEffectsStd{i},anVarsStd{i});
    else
        cfg = buildDesignOneWay(cfg,nSubj,anVarsStd{i});
    end
    cfg.whicheffect = anEffectsStd{i};
    if strcmp(mode,'eeg_sim')
        cfg.numrandomization = 1000;
    else
        cfg.numrandomization = 5000;
    end
    if strcmp(mode,'behav')
        cfg.correctm = 'no';
    else
        cfg.correctm = 'cluster'; %'max'
        cfg.clusterstatistic = 'maxsum';
        cfg.clusterthreshold = 'nonparametric_common';
        cfg.clusteralpha = 0.05;
        cfg.clustertail = 1;
    end
    cfg.alpha = 0.05;
    cfg.tail = 1;
    cfg.randomseed = 'yes';
    if ~isempty(regexp(mode,'sim','once'))
        cfg.feedback = 'no';
    else
        cfg.feedback = 'text';
    end
    cfg.statistic = 'ft_statfun_circ_awtest';
    % Computing statistics 
    tempStats{i} = ft_timelockstatistics(cfg,data{:});
end
if ~strcmp(mode,'behav')
    genTimeStr = 'tr_';
else
    genTimeStr = '';
end

for i = 1:numel(anEffects)
    pName = sprintf('p_%s%s',genTimeStr,anEffects{i});
    hName = sprintf('h_%s%s',genTimeStr,anEffects{i});
    stName = sprintf('st_%s%s',genTimeStr,anEffects{i});
    stats.(pName) = NaN(size(time,2),1);
    stats.(pName)(isInTimeWin) = tempStats{i}.prob';
    stats.(hName) = NaN(size(time,2),1);
    stats.(hName)(isInTimeWin) = tempStats{i}.mask';
    stats.(stName) = NaN(size(time,2),1);
    stats.(stName)(isInTimeWin) = tempStats{i}.stat';
end

end

% Accessory function for building the 3-way design
function cfg = buildDesignThreeWay(cfg,nSubj,anName,varNames)
% varNames must be defined in order A_B_C
nCond = size(varNames,2);
% Split variable names to get individual factor levels
varConds = cellfun(@strsplit,varNames,repmat({'_'},size(varNames)),'UniformOutput',false);
% nConditions x nFactors
varConds = cat(1,varConds{:});
% Independent variables. One row for each factor, recoded to [1,2]
% within each level
ivar = ones(size(varConds));
ivar(ismember(varConds,{'A2','B2','C2'})) = 2;
ivar = repmat(ivar',1,nSubj);
% Unit variable. This indicates the unit of observation, the
% subject. Values are 1:nSubj. This restricts permutations within 
% each subject. 
uvar = repmat(1:nSubj,nCond,1);
uvar = uvar(:)';
% The second dimension restricts the permutations to within certain
% conditions based on the analysis at hand. 
if ismember(anName,{'A','B','C'})
    % Main effects. Permute separetely within the levels of the not
    % tested effects. 
    if strcmp(anName,'A')
        % Test main effect of A, permute within levels of B and C.
        idx = ismember(varConds(:,1),'A1')';
    elseif strcmp(anName,'B')
        % Test main effect of B, permute within levels of A and C.
        idx = ismember(varConds(:,2),'B1')';
    elseif strcmp(anName,'C')
        % Test main effect of C, permute within levels of A and B.
        idx = ismember(varConds(:,3),'C1')';
    end
    % Control variable. The different levels of the control variable 
    % indicate the blocks in which the resampling can be done
    % the replications should not be resampled over the blocks
    temp = NaN(size(idx));
    temp(idx) = 1:sum(idx);
    temp(~idx) = 1:sum(~idx);
    cvar = repmat(temp,1,nSubj);
    % Saving design and related info
    cfg.design = cat(1,uvar,cvar,ivar);
    % cfg.uvar = 1;
    % I treat uvar as cvar
    cfg.cvar = 1:2;
    cfg.ivar = 3:size(cfg.design,1);
elseif ismember(anName,{'A X B','A X C','B X C'}) 
    % Two-way interactions. Permute the simple effects within
    % levels of 3rd factors (Edgington 1995 Randomization tests)
    if strcmp(anName,'A X B')
        % permute along A or B (both equivalent) within
        % levels of C
        cidx = ismember(varConds(:,3),'C1')';
        widx = ismember(varConds(:,2),'B1')';
    elseif strcmp(anName,'A X C')
        % permute along A or C (both equivalent) within
        % levels of B
        cidx = ismember(varConds(:,2),'B1')';
        widx = ismember(varConds(:,3),'C1')';
    elseif strcmp(anName,'B X C')
        % permute along B or C (both equivalent) within
        % levels of A
        cidx = ismember(varConds(:,1),'A1')';
        widx = ismember(varConds(:,2),'B1')';
    end
    % Permutations are restricted within the levels of the control
    % variable
    temp = ones(size(cidx));
    temp(~cidx) = 2;
    cvar = repmat(temp,1,nSubj);
    % Examples within the levels of this variable are kept together
    % during permutations
    temp = ones(size(widx));
    temp(~widx) = 2;
    wvar = repmat(temp,1,nSubj);
    % Saving design and related info
    cfg.design = cat(1,uvar,cvar,wvar,ivar);
    % cfg.uvar = 1;
    % I treat uvar as cvar
    cfg.cvar = 1:2;
    cfg.wvar = 3;
    cfg.ivar = 4:size(cfg.design,1);
elseif strcmp(anName,'A X B X C')
    % Three-way interaction, permute freely, so no further
    % variables need to be defined
    % Saving design and related info
    cfg.design = cat(1,uvar,ivar);
    % I treat uvar as cvar
    cfg.cvar = 1;
    cfg.ivar = 2:size(cfg.design,1);
    
else
    error('mvpa:compAVmodelWeightStats:unknownEffect',...
          'This effect is not implemented!')
end

end

% Accessory function for building the 2-way design
function cfg = buildDesignTwoWay(cfg,nSubj,anName,varNames)
% varNames must be defined in order A_B
nCond = size(varNames,2);
% Split variable names to get individual factor levels
varConds = cellfun(@strsplit,varNames,repmat({'_'},size(varNames)),'UniformOutput',false);
% nConditions x nFactors
varConds = cat(1,varConds{:});
% Independent variables. One row for each factor, recoded to [1,2]
% within each level
ivar = ones(size(varConds));
ivar(ismember(varConds,{'A2','B2'})) = 2;
ivar = repmat(ivar',1,nSubj);
% Unit variable. This indicates the unit of observation, the
% subject. Values are 1:nSubj. This restricts permutations within 
% each subject. 
uvar = repmat(1:nSubj,nCond,1);
uvar = uvar(:)';
% The second dimension restricts the permutations to within certain
% conditions based on the analysis at hand. 
if ismember(anName,{'A','B'})
    % Main effects. Permute separetely within the levels of the not
    % tested effects. 
    if strcmp(anName,'A')
        % Test main effect of A, permute within levels of B. 
        idx = ismember(varConds(:,1),'A1')';
    elseif strcmp(anName,'B')
        % Test main effect of B, permute within levels of A.
        idx = ismember(varConds(:,2),'B1')';
    end
    % Control variable. The different levels of the control variable 
    % indicate the blocks in which the resampling can be done
    % the replications should not be resampled over the blocks
    temp = NaN(size(idx));
    temp(idx) = 1:sum(idx);
    temp(~idx) = 1:sum(~idx);
    cvar = repmat(temp,1,nSubj);
    % Saving design and related info
    cfg.design = cat(1,uvar,cvar,ivar);
    % cfg.uvar = 1;
    % I treat uvar as cvar
    cfg.cvar = 1:2;
    cfg.ivar = 3:size(cfg.design,1);
    
    % A two-way interaction cannot be tested in a 2 X 2 design
    % using this permutation approach. 
    
else
    error('mvpa:compAVmodelWeightStats:unknownEffect',...
          'This effect is not implemented!')
end

end

% Accessory function for building the 3-way design
function cfg = buildDesignOneWay(cfg,nSubj,varNames)
% varNames must contain only A1 and A2
nCond = size(varNames,2);
% Split variable names to get individual factor levels
varConds = cellfun(@strsplit,varNames,repmat({'_'},size(varNames)),'UniformOutput',false);
% nConditions x nFactors
varConds = cat(1,varConds{:});
% Independent variables. One row for each factor, recoded to [1,2]
% within each level
ivar = ones(size(varConds));
ivar(ismember(varConds,{'A2'})) = 2;
ivar = repmat(ivar',1,nSubj);
% Unit variable. This indicates the unit of observation, the
% subject. Values are 1:nSubj. This restricts permutations within 
% each subject. 
uvar = repmat(1:nSubj,nCond,1);
uvar = uvar(:)';
% Since there is only one factor with two levels we can permute
% freely within the unit variable, so no further variables need to
% be defined. 
%Saving design and related info
cfg.design = cat(1,uvar,ivar);
% cfg.uvar = 1;
% I treat uvar as cvar
cfg.cvar = 1;
cfg.ivar = 2;

end
