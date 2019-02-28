function subjects = analysis_behav(varargin)
% Performs detailed behavioural data analysis
%
% USAGE: 
%   subjects = analysis_behav(varargin)
% INPUT: 
%   Optional
%       plotfigures (logical): Whether to plot figures 
% OUTPUT: 
%   subjects (struct): structure arrays of tables containing the recoded 
%       data of individual subjects. 
%

% Copyright(C) 2018, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
checkType = @(x) iscellstr(x) && numel(x) == 1 && all(ismember(x,validTypes));
addOptional(p,'plotfigures',false,@(x) validateattributes(x,...
                                                  {'logical'},{'vector','numel',1}));
parse(p,varargin{:});
plotfigures = p.Results.plotfigures;

%% Collecting data. 
expStage = 'final';
load('subject_spec.mat');
subjList = [subject_spec.subID]';
subjList = cellstr(num2str(subjList));
subjects = struct;

for i = 1:numel(subjList)

    saveDf = cd(DEC_2_setupdir(expStage,'anal_behav_sub',subjList{i}));
    fileList = dir;
    fileList = {fileList.name}';
    matchID = ~cellfun(@isempty,regexp(fileList,'preproc_BEHAV_[0-9]{3}.mat'));
    if sum(matchID) == 0
        cd(saveDf);
        error(['No behavioural file for subject %s! ',...
              'Please run behavioural preprocessing for all subjects'],subjList{i});
    end
    subjects(i).subID = str2double(subjList{i});
    temp = load(fileList{matchID});
    subjects(i).behavData = temp.behavData;
    
    cd(saveDf);
end

smallDisp = 6.67;

%% Computing subject level indices. 
for i = 1:numel(subjects)
    
    % Unisensroy localization correlation
    subjects(i).behavData.corrInput = [subjects(i).behavData.locA subjects(i).behavData.resp]; 
    temp = varfun(@correlatedata,subjects(i).behavData(...
         ~isnan(subjects(i).behavData.resp) & ...    % no NaNs
         isnan(subjects(i).behavData.locV) & ...     % Auditory only trials
         subjects(i).behavData.toReject == 0,:),...  % Valid trials only
         'InputVariables','corrInput');
    subjects(i).audLocCorr = table();
    subjects(i).audLocCorr.rho = temp.correlatedata_corrInput(1);
    subjects(i).audLocCorr.p = temp.correlatedata_corrInput(2);
    subjects(i).audLocCorr.Properties.RowNames = {};
    subjects(i).audLocCorr.subID = repmat(subjects(i).subID,size(subjects(i).audLocCorr.rho));
    subjects(i).audLocCorr = subjects(i).audLocCorr(:,[end,1:end-1]);    
    
    subjects(i).behavData.corrInput = [subjects(i).behavData.locV subjects(i).behavData.resp];
    temp = varfun(@correlatedata,subjects(i).behavData(...
        ~isnan(subjects(i).behavData.resp) & ...    % no NaNs
        isnan(subjects(i).behavData.locA) & ...     % Visual only trials
        subjects(i).behavData.toReject == 0,:),...  % Valid trials only
        'InputVariables','corrInput','GroupingVariables',{'relV'});
    subjects(i).visLocCorr = table();
    subjects(i).visLocCorr.relV = temp.relV;
    subjects(i).visLocCorr.rho = temp.correlatedata_corrInput(:,1);
    subjects(i).visLocCorr.p = temp.correlatedata_corrInput(:,2);
    subjects(i).visLocCorr.Properties.RowNames = {};
    subjects(i).visLocCorr.subID = repmat(subjects(i).subID,size(subjects(i).visLocCorr.rho));
    subjects(i).visLocCorr = subjects(i).visLocCorr(:,[end,1:end-1]);
    subjects(i).behavData.corrInput = [];
    
    % Unisensory localization regression
    subjects(i).behavData.regrInput = [subjects(i).behavData.locA subjects(i).behavData.resp];
    subjects(i).audLocMdl = varfun(@regressdata,subjects(i).behavData(...
                                    ~isnan(subjects(i).behavData.resp) & ...    % no NaNs
                                    isnan(subjects(i).behavData.locV) & ...     % Auditory only trials
                                    subjects(i).behavData.toReject == 0,:),...  % Valid trials only
                                    'InputVariables','regrInput');
    subjects(i).audLocMdl.Properties.VariableNames{'regressdata_regrInput'} = 'mdl';
    subjects(i).audLocMdl.subID = repmat(subjects(i).subID,size(subjects(i).audLocMdl.mdl));
    subjects(i).audLocMdl = subjects(i).audLocMdl(:,[end,1:end-1]);
    subjects(i).audLocMdl.Properties.RowNames = {};
    
    subjects(i).behavData.regrInput = [subjects(i).behavData.locV subjects(i).behavData.resp];
    subjects(i).visLocMdl = varfun(@regressdata,subjects(i).behavData(...
                                    ~isnan(subjects(i).behavData.resp) & ...    % no NaNs
                                    isnan(subjects(i).behavData.locA) & ...     % Visual only trials
                                    subjects(i).behavData.toReject == 0,:),...  % Valid trials only
                                    'InputVariables','regrInput');
    subjects(i).visLocMdl.Properties.VariableNames{'regressdata_regrInput'} = 'mdl';
    subjects(i).visLocMdl.subID = repmat(subjects(i).subID,size(subjects(i).visLocMdl.mdl));
    subjects(i).visLocMdl = subjects(i).visLocMdl(:,[end,1:end-1]);
    subjects(i).visLocMdl.Properties.RowNames = {};
    
    subjects(i).behavData.regrInput = [];
    
    % Unisensroy localization distribution
    locLevels = unique(subjects(i).behavData.locA(~isnan(subjects(i).behavData.locA)))';
    subjects(i).behavData.distribInput = [subjects(i).behavData.resp,repmat(locLevels,size(subjects(i).behavData.resp,1),1)];
    subjects(i).uniLocDistrib = varfun(@calcrespdistrib,subjects(i).behavData(...
        ~isnan(subjects(i).behavData.resp) & ...    % no NaNs
        isnan(subjects(i).behavData.locV) & ...     % Auditory only trials
        subjects(i).behavData.toReject == 0,:),...  % Valid trials only
            'InputVariables','distribInput','GroupingVariables',{'locA'});
    subjects(i).uniLocDistrib.Properties.VariableNames{'calcrespdistrib_distribInput'} = 'distribution';
    subjects(i).uniLocDistrib.Properties.RowNames = {};
    subjects(i).uniLocDistrib.subID = repmat(subjects(i).subID,size(subjects(i).uniLocDistrib.distribution,1),1);
    subjects(i).uniLocDistrib.locV = NaN(size(subjects(i).uniLocDistrib.distribution,1),1);
    subjects(i).uniLocDistrib.relV = NaN(size(subjects(i).uniLocDistrib.distribution,1),1);
    subjects(i).uniLocDistrib = subjects(i).uniLocDistrib(:,[4,1,5,6,2,3]);
    
    temp = varfun(@calcrespdistrib,subjects(i).behavData(...
        ~isnan(subjects(i).behavData.resp) & ...    % no missing responses
        isnan(subjects(i).behavData.locA) & ...     % Visual only trials
        subjects(i).behavData.toReject == 0,:),...  % Valid trials only
            'InputVariables','distribInput','GroupingVariables',{'locV','relV'});
    temp.Properties.VariableNames{'calcrespdistrib_distribInput'} = 'distribution';
    temp.Properties.RowNames = {};
    temp.subID = repmat(subjects(i).subID,size(temp.distribution,1),1);
    temp.locA = NaN(size(temp.distribution,1),1);
    temp = temp(:,[5:6,1:4]);
    temp = sortrows(temp,{'relV','locV'});
    
    subjects(i).uniLocDistrib = [subjects(i).uniLocDistrib;temp];
    
    % Bisensory localization distribution
    subjects(i).biLocDistrib = varfun(@calcrespdistrib,subjects(i).behavData(...
        ~isnan(subjects(i).behavData.resp) & ...    % no NaNs
        ~isnan(subjects(i).behavData.locV) & ...    % Audio-visual trials
        ~isnan(subjects(i).behavData.locA) & ...
        subjects(i).behavData.toReject == 0,:),...  % Valid trials only
            'InputVariables','distribInput','GroupingVariables',...
            {'locA','locV','relV','task'});
    subjects(i).biLocDistrib.Properties.VariableNames{'calcrespdistrib_distribInput'} = 'distribution';
    subjects(i).biLocDistrib.Properties.RowNames = {};
    subjects(i).biLocDistrib.subID = repmat(subjects(i).subID,size(subjects(i).biLocDistrib.distribution,1),1);
    subjects(i).biLocDistrib = subjects(i).biLocDistrib(:,[7,1:6]);
    
    subjects(i).behavData.distribInput = [];
    
    % Calculating absolute visual bias
    subjects(i).behavData.disparity = subjects(i).behavData.locV - subjects(i).behavData.locA;
    subjects(i).behavData.disparity(subjects(i).behavData.disparity > 6.65 & ...
                                 subjects(i).behavData.disparity < 6.68) = smallDisp;
    subjects(i).behavData.disparity(subjects(i).behavData.disparity < -6.65 & ...
                                 subjects(i).behavData.disparity > -6.68) = -smallDisp;
    subjects(i).behavData.absDisparity = abs(subjects(i).behavData.disparity);
    subjects(i).behavData.absVisBias = subjects(i).behavData.resp - subjects(i).behavData.locA;
    subjects(i).behavData.absVisBias(subjects(i).behavData.task == 2) = NaN;    % just for auditory localization
    subjects(i).behavData.absVisBias(isnan(subjects(i).behavData.locV)) = NaN;  % just for AV trials
    subjects(i).absVisBias = join(varfun(@mean,subjects(i).behavData(...
        ~isnan(subjects(i).behavData.resp) & ...     % no missing responses
        subjects(i).behavData.toReject == 0 & ...    % Valid trials only
        ~isnan(subjects(i).behavData.locA) & ...    % AV trials
        ~isnan(subjects(i).behavData.locV) & ...
        subjects(i).behavData.task == 1,:),...      % auditory localization
        'InputVariables','absVisBias','GroupingVariables',{'relV','disparity'}),...
        varfun(@std,subjects(i).behavData(...
        ~isnan(subjects(i).behavData.resp) & ...     % no missing responses
        subjects(i).behavData.toReject == 0 & ...    % Valid trials only
        ~isnan(subjects(i).behavData.locA) & ...    % AV trials
        ~isnan(subjects(i).behavData.locV) & ...
        subjects(i).behavData.task == 1,:),...      % auditory localization
        'InputVariables','absVisBias','GroupingVariables',{'relV','disparity'}));
    subjects(i).absVisBias.Properties.VariableNames{'mean_absVisBias'} = 'mean';
    subjects(i).absVisBias.Properties.VariableNames{'std_absVisBias'} = 'std';
    subjects(i).absVisBias.sem = subjects(i).absVisBias.std./sqrt(subjects(i).absVisBias.GroupCount);
    subjects(i).absVisBias.subID = repmat(subjects(i).subID,size(subjects(i).absVisBias.mean));
    subjects(i).absVisBias = subjects(i).absVisBias(:,[7 1:6]);
    subjects(i).absVisBias.Properties.RowNames = {};
    
    % Modeling responses
    groupings = {'all','cmb_low','cmb_high','psalpha_low','psalpha_high'};
    betas = table;
    for j = 1:size(groupings,2)
        % Congruent trials are ignored to be consistent with the EEG analysis
        subjects(i).behavData.largeAVdisp = subjects(i).behavData.absDisparity > smallDisp;
        subjects(i).behavData.smallAVdisp = subjects(i).behavData.absDisparity <= smallDisp;

        if strcmp(groupings{j},'all')
            isExample = ~isnan(subjects(i).behavData.resp) & ... % no missing responses
                subjects(i).behavData.toReject == 0 & ...         % Valid trials onlys
                ~isnan(subjects(i).behavData.locA) & ...          % AV trials
                ~isnan(subjects(i).behavData.locV);
            
            locGLM = fitparams(subjects(i).behavData(isExample,:));
            % Saving betas in a convenient fashion
            temp = locGLM.Coefficients.Estimate(2:end)';
            temp = array2table(temp,'VariableNames',locGLM.Coefficients.Properties.RowNames(2:end));
            temp.subID = subjects(i).subID;
            temp.groupings = groupings(j);
            betas = cat(1,betas,temp);
            
        elseif regexp(groupings{j},'^cmb_','once')
            isExample = ~isnan(subjects(i).behavData.resp) & ... % no missing responses
                subjects(i).behavData.toReject == 0 & ...         % Valid trials onlys
                ~isnan(subjects(i).behavData.locA) & ...          % AV trials
                ~isnan(subjects(i).behavData.locV) & ...
                ~isnan(subjects(i).behavData.cmb);
            if strcmp(groupings{j},'cmb_low')
                isExample = isExample & subjects(i).behavData.cmb <= 0;
            else
                isExample = isExample & subjects(i).behavData.cmb > 0;
            end
            
            locGLM = fitparams(subjects(i).behavData(isExample,:));
            % Saving betas in a convenient fashion
            temp = locGLM.Coefficients.Estimate(2:end)';
            temp = array2table(temp,'VariableNames',locGLM.Coefficients.Properties.RowNames(2:end));
            temp.subID = subjects(i).subID;
            temp.groupings = groupings(j);
            betas = cat(1,betas,temp);
            
        elseif regexp(groupings{j},'^psalpha_','once')
            
            psalphaTime = subjects(i).behavData.Properties.UserData.psalpha_cols_time;
            timeWin = {[-0.6,-0.5],[-0.5,-0.4],[-0.4,-0.3],[-0.3,-0.2],[-0.2,-0.1],...
                      [-0.1,-0]};
            fname = 'psalpha_bin_timewins';
            if ~isfield(subjects(i).behavData.Properties.UserData,fname)
                subjects(i).behavData.Properties.UserData.(fname) = ...
                    timeWin;
            end
            
            for iTimeWin = 1:numel(timeWin)
                isExample = ~isnan(subjects(i).behavData.resp) & ... % no missing responses
                subjects(i).behavData.toReject == 0 & ...         % Valid trials onlys
                ~isnan(subjects(i).behavData.locA) & ...          % AV trials
                ~isnan(subjects(i).behavData.locV);
                
                % Find time window of interest
                isTimeWin = psalphaTime >= timeWin{iTimeWin}(1) & ...
                    psalphaTime < timeWin{iTimeWin}(2);
                % Calculate mean alpha power over the timewindow
                subjects(i).behavData.psalpha_temp = NaN(size(subjects(i).behavData,1),1);
                subjects(i).behavData.psalpha_temp(isExample) = ...
                    nanmean(subjects(i).behavData.psalpha(isExample,isTimeWin),2);
                % Calculate median alpha power separately for each condition
                m = varfun(@nanmedian,subjects(i).behavData(isExample,:),...
                           'InputVariables',{'psalpha_temp'},...
                           'GroupingVariables',{'condition'});
                temp = NaN(size(subjects(i).behavData,1),1);
                for cond = m.condition'
                    mCond = m.nanmedian_psalpha_temp(m.condition == cond);
                    temp(subjects(i).behavData.condition == cond & ...
                         subjects(i).behavData.psalpha_temp <= mCond) = 0;
                    temp(subjects(i).behavData.condition == cond & ...
                         subjects(i).behavData.psalpha_temp > mCond) = 1;
                end
            
                subjects(i).behavData.psalpha_bin(:,iTimeWin) = temp;
                if strcmp(groupings{j},'psalpha_low')
                    isExample = isExample & subjects(i).behavData.psalpha_bin(:,iTimeWin) == 0;
                else
                    isExample = isExample & subjects(i).behavData.psalpha_bin(:,iTimeWin) == 1;
                end
                
                locGLM = fitparams(subjects(i).behavData(isExample,:));
                % Saving betas in a convenient fashion
                temp = locGLM.Coefficients.Estimate(2:end)';
                temp = array2table(temp,'VariableNames',locGLM.Coefficients.Properties.RowNames(2:end));
                temp.subID = subjects(i).subID;
                temp.groupings = {[groupings{j},'_',num2str(iTimeWin)]};
                betas = cat(1,betas,temp);
            end
        end
        
    end
    betas.groupings = categorical(betas.groupings);
    subjects(i).avBetas = betas;
    nVar = size(subjects(i).avBetas,2);
    subjects(i).avBetas = subjects(i).avBetas(:,[nVar-1,nVar,1:nVar-2]);
    % Computing individual AV weights
    funb2w = @(x) atan2(x(:,1),x(:,2));
    varNames = sort(subjects(i).avBetas.Properties.VariableNames);
    a_varNames = varNames(~cellfun(@isempty,regexp(varNames,'^A_')));
    v_varNames = varNames(~cellfun(@isempty,regexp(varNames,'^V_')));
    wav_varNames = strrep(a_varNames,'A_','');
    % Pairing A and V betas for conditions respectivel.
    temp = table;
    for k = 1:size(wav_varNames,2)
        temp.(wav_varNames{k}) = [subjects(i).avBetas.(v_varNames{k}),subjects(i).avBetas.(a_varNames{k})];
    end
    temp.groupings = subjects(i).avBetas.groupings;
    subjects(i).wav = varfun(funb2w,temp,'GroupingVariables','groupings');
    varNames = subjects(i).wav.Properties.VariableNames;
    isWavVar = varNames(~cellfun(@isempty,regexp(varNames,'^Fun_')));
    subjects(i).wav.Properties.VariableNames(isWavVar) = strcat(wav_varNames,repmat({'_wav'},size(wav_varNames)));
    subjects(i).wav.subID = repmat(subjects(i).subID,size(subjects(i).wav,1),1);
    subjects(i).wav.GroupCount = [];
    subjects(i).wav.Properties.RowNames = {};
    nVar = size(subjects(i).wav,2);
    subjects(i).wav = subjects(i).wav(:,[nVar,1:nVar-1]);
    
end

if plotfigures(1)
    % Plotting subject level data.
    for i = 1:numel(subjects)
        
        hFig = figure('Units','Normalized','Position',[0.125,0,0.75,1]);
        % Plotting unisensory localization data
        plotlocdistributions(subjects(i).uniLocDistrib,hFig,subplot(2,2,1,'OuterPosition',[0,0.5,0.55,0.45]));
        
        % Plotting absolute bias vs av-disparity
        % Absolute bias is computed only for AV-a trials.
        plotabsbiasvsdisp(subjects(i).absVisBias,subjects(i).behavData(~isnan(subjects(i).behavData.absVisBias),:),...
            hFig,subplot(2,2,2,'OuterPosition',[0.666,0.5,0.31,0.45]));
        
        % Plotting response times
        plotresptimes(subjects(i).behavData,hFig,subplot(2,2,3,'OuterPosition',[0,0,0.3,0.45]));
        
        % Plotting weights
%         plotweights(subjects(i).locGLM,hFig,subplot(2,2,4,'OuterPosition',[0.333,0,0.65,0.45]));
        
        suplabel(sprintf('Subject %d behavioural results',subjects(i).subID),'t',[0,0,0.95,0.95]);
    end
end
%% Computing group level results

%               Group level bisensory localization distribution
%--------------------------------------------------------------------------
grBiLocDistrib.data = cat(1,subjects.biLocDistrib);
fieldName = 'meanDistribution';
grBiLocDistrib.(fieldName) = varfun(@mean,grBiLocDistrib.data,...
    'InputVariables','distribution','GroupingVariables',...
    {'locA','locV','relV','task'});
grBiLocDistrib.(fieldName).Properties.RowNames = {};
grBiLocDistrib.(fieldName).Properties.VariableNames{'mean_distribution'} = 'distribution';

%                          Group audio visual betas
%--------------------------------------------------------------------------
grBetas.data = cat(1,subjects.avBetas);
varNames = grBetas.data.Properties.VariableNames;
betaVarNames = varNames(~cellfun(@isempty,regexp(varNames,'^(A|V)_')));
fieldName = 'meanBeta';

% Computing mean and standard deviation
ci_alpha = 0.32;
wrap_confmean = @(x) mean(x) - mvpa.confmean(x,1,ci_alpha);
% Mean betas
grBetas.(fieldName) = varfun(@mean,grBetas.data,...
    'InputVariables',betaVarNames,'GroupingVariables','groupings');
grBetas.(fieldName).Properties.RowNames = {};
varNames = grBetas.(fieldName).Properties.VariableNames;
funIdx = ~cellfun(@isempty,regexp(varNames,'^mean_'));
funNames = varNames(funIdx);
grBetas.(fieldName).Properties.VariableNames(funIdx) = strrep(funNames,'mean_','');
% Computing confidence intervals
temp = varfun(wrap_confmean,grBetas.data,...
    'InputVariables',betaVarNames,'GroupingVariables','groupings');
temp.Properties.RowNames = {};
varNames = temp.Properties.VariableNames;
funIdx = ~cellfun(@isempty,regexp(varNames,'^Fun_'));
funNames = varNames(funIdx);
temp.Properties.VariableNames(funIdx) = strcat(strrep(funNames,'Fun_',''),'_ci');
% Joining the tables
grBetas.(fieldName) = join(grBetas.(fieldName),temp);

%                          Group audio visual weights
%--------------------------------------------------------------------------
grWav.data = cat(1,subjects.wav);
varNames = grWav.data.Properties.VariableNames;
wavVarNames = varNames(~cellfun(@isempty,regexp(varNames,'_wav$')));
fieldName = 'meanWav';

% Computing mean and standard deviation
wrap_circ_mean = @(x) degrees(circ_mean(radians(x)));
wrap_circ_confmean = @(x) degrees(circ_confmean(radians(x),ci_alpha));
% Computing circular means
grWav.(fieldName) = varfun(wrap_circ_mean,grWav.data,...
    'InputVariables',wavVarNames,'GroupingVariables','groupings');
grWav.(fieldName).Properties.RowNames = {};
varNames = grWav.(fieldName).Properties.VariableNames;
funIdx = ~cellfun(@isempty,regexp(varNames,'^Fun_'));
funNames = varNames(funIdx);
grWav.(fieldName).Properties.VariableNames(funIdx) = strrep(funNames,'Fun_','');
% Computing confidence intervals
temp = varfun(wrap_circ_confmean,grWav.data,...
    'InputVariables',wavVarNames,'GroupingVariables','groupings');
temp.Properties.RowNames = {};
varNames = temp.Properties.VariableNames;
funIdx = ~cellfun(@isempty,regexp(varNames,'^Fun_'));
funNames = varNames(funIdx);
temp.Properties.VariableNames(funIdx) = strrep(strrep(funNames,'Fun_',''),'_wav','_ci');
% Joining the tables
grWav.(fieldName) = join(grWav.(fieldName),temp);


%             Group audio and visual localization performance 
%--------------------------------------------------------------------------
% Group regressions
temp = cat(1,subjects.visLocMdl);
temp = cellfun(@(x) x.Coefficients,temp.mdl,'UniformOutput',false);
for i = 1:numel(temp)
    temp{i}.EstimateName = categorical(temp{i}.Properties.RowNames);
    temp{i}.Properties.RowNames = {};
    temp{i}.subID = repmat(subjects(i).subID,size(temp{i},1),1);
    temp{i} = temp{i}(:,[end,end-1,1:end-2]);
end
grVisLocMdl.data = cat(1,temp{:});
% Computing mean and standard deviation
fieldName = 'meanEstimates';
grVisLocMdl.(fieldName) = join(varfun(@mean,grVisLocMdl.data,...
    'InputVariables','Estimate','GroupingVariables','EstimateName'),...
    varfun(@std,grVisLocMdl.data,...
    'InputVariables','Estimate','GroupingVariables','EstimateName'));
grVisLocMdl.(fieldName).Properties.VariableNames{'mean_Estimate'} = 'mean';
grVisLocMdl.(fieldName).Properties.VariableNames{'std_Estimate'} = 'std';
grVisLocMdl.(fieldName).sem = grVisLocMdl.(fieldName).std ./ ...
    sqrt(grVisLocMdl.(fieldName).GroupCount);

temp = cat(1,subjects.audLocMdl);
temp = cellfun(@(x) x.Coefficients,temp.mdl,'UniformOutput',false);
for i = 1:numel(temp)
    temp{i}.EstimateName = categorical(temp{i}.Properties.RowNames);
    temp{i}.Properties.RowNames = {};
    temp{i}.subID = repmat(subjects(i).subID,size(temp{i},1),1);
    temp{i} = temp{i}(:,[end,end-1,1:end-2]);
end
grAudLocMdl.data = cat(1,temp{:});
% Computing mean and standard deviation
fieldName = 'meanEstimates';
grAudLocMdl.(fieldName) = join(varfun(@mean,grAudLocMdl.data,...
    'InputVariables','Estimate','GroupingVariables','EstimateName'),...
    varfun(@std,grAudLocMdl.data,...
    'InputVariables','Estimate','GroupingVariables','EstimateName'));
grAudLocMdl.(fieldName).Properties.VariableNames{'mean_Estimate'} = 'mean';
grAudLocMdl.(fieldName).Properties.VariableNames{'std_Estimate'} = 'std';
grAudLocMdl.(fieldName).sem = grAudLocMdl.(fieldName).std ./ ...
    sqrt(grAudLocMdl.(fieldName).GroupCount);

% Group correlations
grVisLocCorr.data = cat(1,subjects.visLocCorr);
% Checking if any of the subject have rho = 1, and changing them to 0.9999
% to avoid Inf later on. 
grVisLocCorr.data.rho(grVisLocCorr.data.rho == 1) = 0.9999;
% Using a dummy variable to leverage grouping feature of varfun
% grVisLocCorr.dxata.dummy = ones(size(grVisLocCorr.data,1),1);
% Computing mean and standard deviation after Fisher transform,
% then inverse Fisher transform
meanFisherTrFun = @(x) tanh(nanmean(atanh(x)));
stdFisherTrFun = @(x) tanh(nanstd(atanh(x)));
fieldName = 'meanCorrelation';
grVisLocCorr.(fieldName) = varfun(meanFisherTrFun,...
    grVisLocCorr.data,'InputVariables',{'rho'},'GroupingVariables',{'relV'});
grVisLocCorr.(fieldName).Properties.VariableNames{'Fun_rho'} = 'mean';
grVisLocCorr.(fieldName) = join(grVisLocCorr.(fieldName),...
    varfun(stdFisherTrFun,grVisLocCorr.data,...
    'InputVariables',{'rho'},'GroupingVariables',{'relV'}));
grVisLocCorr.(fieldName).Properties.VariableNames{'Fun_rho'} = 'std';
grVisLocCorr.(fieldName).sem = grVisLocCorr.(fieldName).std ./ ...
    sqrt(grVisLocCorr.(fieldName).GroupCount);
% grVisLocCorr.(fieldName).dummy = [];
% grVisLocCorr.data.dummy = [];

grAudLocCorr.data = cat(1,subjects.audLocCorr);
% Using a dummy variable to leverage grouping feature of varfun
grAudLocCorr.data.dummy = ones(size(grAudLocCorr.data,1),1);
% Computing mean and standard deviation after Fisher transform,
% then inverse Fisher transform
fieldName = 'meanCorrelation';
grAudLocCorr.(fieldName) = varfun(meanFisherTrFun,...
    grAudLocCorr.data,'InputVariables','rho','GroupingVariables','dummy');
grAudLocCorr.(fieldName).Properties.VariableNames{'Fun_rho'} = 'mean';
grAudLocCorr.(fieldName) = join(grAudLocCorr.(fieldName),...
    varfun(stdFisherTrFun,grAudLocCorr.data,...
    'InputVariables','rho','GroupingVariables','dummy'));
grAudLocCorr.(fieldName).Properties.VariableNames{'Fun_rho'} = 'std';
grAudLocCorr.(fieldName).sem = grAudLocCorr.(fieldName).std ./ ...
    sqrt(grAudLocCorr.(fieldName).GroupCount);
grAudLocCorr.(fieldName).dummy = [];
grAudLocCorr.data.dummy = [];

% clear badLocList i mdl xData yData temp toBeExcluded

%% Saving individual behavioural results
fprintf('\n\nSaving data...\n\n');
for i = 1:size(subjects,2)
    dataBehav = subjects(i);
    subID = num2str(dataBehav.subID);
    savePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),['BEHAV_ANAL_',subID,'.mat']);
    save(savePath,'-struct','dataBehav','-v7.3');
end

savePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_group'),'BEHAV_ANAL_group.mat');
if exist(savePath,'file')
    warning(['Group level behavioural analysis file already exists',...
             ', it will be overwritten']);
end

save(savePath,'subjects','grAudLocMdl','grVisLocMdl', ...
     'grAudLocCorr','grVisLocCorr','grBetas','grWav','-v7.3');

end

%% Auxilliary functions
function out = correlatedata(X)

[rho,p] = corr(X);

out = [rho(2,1),p(2,1)];

end


function out = regressdata(X)

mdl = fitglm(X(:,1:end-1),X(:,end),'linear','Distribution','normal');

out = {mdl};
end


function out = calcrespdistrib(X)

resp = X(:,1);
levels = X(1,2:end);
nLevels = numel(levels);
out = zeros(1,nLevels);

for i = 1:nLevels
    out(i) = sum(resp == levels(i))/size(resp,1);
end

end


function mdl = fitparams(behavData)
% Summary of this function goes here
% Detailed explanation goes here

respLabels = behavData.resp;
nExamples = size(behavData,1);

r_d_a = behavData.relV == 12 & behavData.smallAVdisp & behavData.task == 1;
r_d_v = behavData.relV == 12 & behavData.smallAVdisp & behavData.task == 2;
r_D_a = behavData.relV == 12 & behavData.largeAVdisp & behavData.task == 1;
r_D_v = behavData.relV == 12 & behavData.largeAVdisp & behavData.task == 2;
R_d_a = behavData.relV == 2 & behavData.smallAVdisp & behavData.task == 1;
R_d_v = behavData.relV == 2 & behavData.smallAVdisp & behavData.task == 2;
R_D_a = behavData.relV == 2 & behavData.largeAVdisp & behavData.task == 1;
R_D_v = behavData.relV == 2 & behavData.largeAVdisp & behavData.task == 2;

variables = zeros(nExamples,16);
variableNames = cell(1,17);
variables(r_d_a,1) = behavData.locA(r_d_a); variableNames{1} = 'A_r_d_a';
variables(r_d_v,2) = behavData.locA(r_d_v); variableNames{2} = 'A_r_d_v';
variables(r_D_a,3) = behavData.locA(r_D_a); variableNames{3} = 'A_r_D_a';
variables(r_D_v,4) = behavData.locA(r_D_v); variableNames{4} = 'A_r_D_v';
variables(r_d_a,5) = behavData.locV(r_d_a); variableNames{5} = 'V_r_d_a';
variables(r_d_v,6) = behavData.locV(r_d_v); variableNames{6} = 'V_r_d_v';
variables(r_D_a,7) = behavData.locV(r_D_a); variableNames{7} = 'V_r_D_a';
variables(r_D_v,8) = behavData.locV(r_D_v); variableNames{8} = 'V_r_D_v';
variables(R_d_a,9) = behavData.locA(R_d_a); variableNames{9} = 'A_R_d_a';
variables(R_d_v,10) = behavData.locA(R_d_v); variableNames{10} = 'A_R_d_v';
variables(R_D_a,11) = behavData.locA(R_D_a); variableNames{11} = 'A_R_D_a';
variables(R_D_v,12) = behavData.locA(R_D_v); variableNames{12} = 'A_R_D_v';
variables(R_d_a,13) = behavData.locV(R_d_a); variableNames{13} = 'V_R_d_a';
variables(R_d_v,14) = behavData.locV(R_d_v); variableNames{14} = 'V_R_d_v';
variables(R_D_a,15) = behavData.locV(R_D_a); variableNames{15} = 'V_R_D_a';
variables(R_D_v,16) = behavData.locV(R_D_v); variableNames{16} = 'V_R_D_v';
variableNames{17} = 'resp';
mdl = fitglm(variables,respLabels,'linear','Distribution','normal',...
    'VarNames',variableNames);


end

