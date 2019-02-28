function [feat,info,misc] = extractdata(dataDir,condDef,idStr,tr_method)
% Extracts the necessary data for MVPA for a given subject. 
%
% Input: 
%   dataDir: the full path of the data folder to work on. 
%   condDef: table which contains the necessary information to define the
%       conditions. 
%       Variables (in order): 
%           condition task locationAuditory locationVisual
%   idStr: The identification string which determines which eeg files are
%       to be used. This is the first part of the name of the file ahead of
%       the first '_' character. 
%   
% Output: 
%   info: table containing all information necessary from the examples
%       for MVPA. The size is the number of examples x 5 variables
%       Variables (in order): 
%           session
%           condition
%           task
%           locA
%           locV
%   feat: L x M x N 3D or L x M x N x O 4D array containing the features, 
%       if 3D: 
%           L is the number of channels(features)
%           M is the number of time samples
%           N is the number of examples 
%       if 4D: 
%           L is the number of channels(features)
%           M is the number of frequencies
%           N is the number of time samples
%           O is the number of examples 
%
% Details: 
%   Trials marked as 'bad' during the preprocessing are ignored. 
%   size(data.info,1) == N should be true! 
%
%% Parsing input
p = inputParser;

validTrMethods = {'sample-wise','sample-wise-bp','sample-wise-rl',...
                  'sample-wise-sm','sample-wise-source','sample-wise-tf'};

addRequired(p,'dataDir',@ischar);
addRequired(p,'condDef',@istable);
addRequired(p,'idStr',@ischar);
addRequired(p,'tr_method',@(x)any(validatestring(x,validTrMethods)));

parse(p,dataDir,condDef,idStr,tr_method);

dataDir = p.Results.dataDir;
condDef = p.Results.condDef;
idStr = p.Results.idStr;
tr_method = p.Results.tr_method;

%%
saveDF = cd(dataDir);

% list matching files
fileNameList = dir;
fileNameList = {fileNameList.name}';
matchID = ~cellfun(@isempty,regexp(fileNameList,idStr));
if sum(matchID) < 1
    error('mvpa:extractdata:invalidInput',...
        'No eeg files were found for feature extraction! ');
else
    fileNameList = fileNameList(matchID);
end

filePathList = fullfile(repmat({pwd},size(fileNameList)),fileNameList);

cd(saveDF);

[infos,feats] = deal(cell(1,size(fileNameList,1)));

for i = 1:size(filePathList,1)
    
    if ~strcmp(tr_method,'sample-wise-source')
        % Load fieldtrip data file
        ftData = load(filePathList{i});
        fieldNames = fieldnames(ftData);
        if numel(fieldNames) == 1
            ftData = ftData.(fieldNames{1});
        else
            error('mvpa:extractdata:invalidInput',...
                  'Data files must a single fieldtrip data structure');
        end
        
        % get trial info
        if ~isfield(ftData,'trialinfo')
            error('mvpa:extractdata:missingField',...
                  'FieldTrip data must contain the trialinfo field');
        end
        trialInfo = ftData.trialinfo;
        
        % logical vector for good trials
        isGoodTrial = true(size(trialInfo,1),1);
        isGoodTrial(trialInfo(:,end) ~= 0) = false;
        
        % getting the features
        if strcmp(ft_datatype(ftData),'raw')
            tempFeats = cat(3,ftData.trial{:});
            tempFeats = tempFeats(:,:,isGoodTrial);
        else
            tempFeats = ftData.powspctrm(:,:,:,isGoodTrial);
        end
    else
        % Load fieldtrip data file
        D = spm_eeg_load(filePathList{i});
        
        % get trial info
        if ~isfield(D,'trialinfo')
            error('mvpa:extractdata:missingField',...
                  'SPM data must contain the trialinfo field');
        end
        trialInfo = D.trialinfo;
        
        % logical vector for good trials
        isGoodTrial = true(size(trialInfo,1),1);
        isGoodTrial(trialInfo(:,end) ~= 0) = false;
        
        % getting the features
        if size(D.size,2) > 3
            tempFeats = D.selectdata(D.chanlabels,[],[],[]);
            tempFeats = tempFeats(:,:,:,isGoodTrial);
        else
            tempFeats = D.selectdata(D.chanlabels,[],[]);
            tempFeats = tempFeats(:,:,isGoodTrial);
        end
    end
    
    % preallocating vectors and reading data from the data structure
    [task,locA,locV,relV] = deal(NaN(size(trialInfo,1),1));
    
    % assigning columns of trialInfo to the appropriate variables
    session = trialInfo(:,1);
    iTrialSession = trialInfo(:,2);
    condition = trialInfo(:,3);
    hand = trialInfo(:,4);
    resp = trialInfo(:,5);
    
    % getting the infos 
    tempInfo = table(session,iTrialSession,condition,hand,task,locA,locV,relV,resp);
    tempInfo = tempInfo(isGoodTrial,:);
    
    % filling up the task, locA, locV and relV variables according to the
    % conditions
    for j = 1:size(condDef.condition,1)
        tempInfo(tempInfo.condition == condDef.condition(j),5:8) = ...
            repmat(condDef(condDef.condition(j),2:5),...
            size(tempInfo.condition(tempInfo.condition == condDef.condition(j)),1),1);
    end
    
    % Adding the AV disparity variable (derivative of lovA and
    % locV) 1: small, 2: large
    tempInfo.disp = ones(size(tempInfo,1),1);
    tempInfo.disp(abs(tempInfo.locA-tempInfo.locV) > mvpa.smallDisp) = 2;
    
    infos{i} = tempInfo;
    feats{i} = tempFeats;
    
    % Free up memory
    tempInfo = [];
    tempFeats = [];
end

if ~strcmp(tr_method,'sample-wise-source')
    if isfield(ftData,'fsample')
        misc.fs = ftData.fsample;
    elseif isfield(ftData,'time')
        misc.fs = round(1/(ftData.time{1,1}(2)-ftData.time{1,1}(2)));
    else
        error('mvpa:extractdata:missingField',...
              'FieldTrip data must contain the time field');
    end

    misc.timeOnset = ftData.time{1,1}(1);
    misc.featureLabels = ftData.label;
else
    misc.fs = D.fsample;
    misc.timeOnset = D.timeonset;
    misc.featureLabels = D.chanlabels;    
end

if strcmp(tr_method,'sample-wise-tf')
    error('mvpa:extractdata:missingFunctionality',...
        'This part has yet to be implemented');
%     history = D.history;
%     tfInd = ismember({history.fun},'spm_eeg_tf');
%     if any(tfInd)
%         misc.freqres = history(tfInd).args.settings.freqres;
%         misc.frequencies = D.frequencies;
%     end
elseif strcmp(tr_method,'sample-wise-bp')
    cfg = ftData.cfg;
    while 1
        if isfield(cfg,'bpfreq')
            misc.frequencies = cfg.bpfreq;
            break;
        elseif isfield(cfg,'previous')
            cfg = cfg.previous;
        else
            error('mvpa:extractdata:missingField',...
                'FieldTrip data must contain the cfg.bpfreq field.');
        end
    end
end
featSize = size(feats{1,1});

info = cat(1,infos{:});
infos = [];

if size(featSize,2) > 3
    feat = cat(4,feats{:});
    feats = [];
    if size(info,1) ~= size(feat,4)
        error('Inconsistent number of examples! ');
    end
else
    feat = cat(3,feats{:});
    feats = [];
    if size(info,1) ~= size(feat,3)
        error('Inconsistent number of examples! ');
    end
end

end
