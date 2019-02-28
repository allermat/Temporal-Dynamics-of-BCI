function analysis_behav_descriptives(subID)
% Parsing input, checking matlab
p = inputParser;
addRequired(p,'subID',@ischar);
parse(p,subID);
subID = p.Results.subID;

% Checking if there is existing processed behavioural data
expStage = 'final';
behavFilePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),...
    ['preproc_BEHAV_',subID,'.mat']);
if exist(behavFilePath,'file')
    behavData = load(behavFilePath);
    behavData = behavData.dataBehav;
else
    warning('analysis_behav_descriptives:missingData',...
        'Missing behavioural data, returning');
    return;
end

% Load condition definition file
condDef = load('condDef.mat');
condDef = condDef.condDef;

% Location levels
locLevels = unique(behavData.locA(~isnan(behavData.locA)))';
% Function for calculating response distribution
    function out = calcrespdistrib(resp)
        nLevels = numel(locLevels);
        out = zeros(1,nLevels);
        for i = 1:nLevels
            out(i) = sum(resp == locLevels(i))/size(resp,1);
        end
    end

% Bisensory localization distribution
respDistribution = varfun(@calcrespdistrib,behavData(...
    ~isnan(behavData.resp) & ...    % no NaNs
    ~isnan(behavData.locV) & ...    % Audio-visual trials
    ~isnan(behavData.locA) & ...
    behavData.toReject == 0,:),...  % Valid trials only
    'InputVariables','resp','GroupingVariables',...
    {'condition'});
respDistribution.Properties.VariableNames{'Fun_resp'} = 'distribution';
respDistribution.Properties.RowNames = {};
% Unisensroy auditory localization distribution
temp = varfun(@calcrespdistrib,behavData(...
    ~isnan(behavData.resp) & ...    % no NaNs
    isnan(behavData.locV) & ...     % Auditory only trials
    behavData.toReject == 0,:),...  % Valid trials only
    'InputVariables','resp','GroupingVariables',{'condition'});
temp.Properties.VariableNames{'Fun_resp'} = 'distribution';
temp.Properties.RowNames = {};
respDistribution = [respDistribution;temp];
% Unisensroy visual localization distribution
temp = varfun(@calcrespdistrib,behavData(...
    ~isnan(behavData.resp) & ...    % no missing responses
    isnan(behavData.locA) & ...     % Visual only trials
    behavData.toReject == 0,:),...  % Valid trials only
    'InputVariables','resp','GroupingVariables',{'condition'});
temp.Properties.VariableNames{'Fun_resp'} = 'distribution';
temp.Properties.RowNames = {};
respDistribution = [respDistribution;temp];
% Making sure everything is sorted according to conditions
respDistribution = sortrows(respDistribution,'condition');
condDef = sortrows(condDef,'condition');
respDistribution = horzcat(respDistribution,...
    condDef(:,{'locationAuditory','locationVisual','reliabilityVisual','task'}));
respDistribution.Properties.VariableNames{'locationAuditory'} = 'locA';
respDistribution.Properties.VariableNames{'locationVisual'} = 'locV';
respDistribution.Properties.VariableNames{'reliabilityVisual'} = 'relV';
respDistribution.subID = repmat(subID,size(respDistribution.distribution,1),1);
respDistribution = respDistribution(:,{'subID','condition','locA','locV',...
    'relV','task','GroupCount','distribution'});

% Saving subject specific behavioural data
fprintf('\n\nSaving data...\n\n');
savePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),...
    ['descriptives_BEHAV_',subID,'.mat']);
save(savePath,'respDistribution','-v7.3');

end

