function subjects = preproc_behav(varargin)
% Collects all the data necessary for evaluation.
%
% USAGE: 
%   subjects = preproc_behav(varargin)
% INPUT: 
%   Optional
%       type (string): data type. Possible values: 'exp','train', default: 
%           'exp' 
% OUTPUT: 
%   subjects (struct): structure arrays of tables containing the recoded 
%       data of individual subjects. 
%

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validTypes = {'exp','train'};
checkType = @(x) iscellstr(x) && numel(x) == 1 && all(ismember(x,validTypes));
addOptional(p,'type',{'exp'},checkType);
parse(p,varargin{:});
type = p.Results.type;

% Defining the data folder of the experiment's data. 
expStage = 'final';
dirBehavData = DEC_2_setupdir(expStage,'data_behav');
% List all subfolders in the data root folder (these should be the
% individual subjects)
listing = dir(dirBehavData);
isFold = [listing(:).isdir];
nameFolds = {listing(isFold).name}';
nameFolds(ismember(nameFolds,{'.','..','excluded'})) = [];

subjects = struct;
index = 1;
for i = 1:length(nameFolds)
    
    subID = nameFolds{i};
    
    % Collecting basic behavioural data
    behavData = collectruns(fullfile(dirBehavData,subID),type{1});
    
    % Skip excluded subjects. 
    if isempty(behavData), continue; end
    
    % Computing crossmodal bias 
    % Selecting incongruent audiovisual trials with auditory response
    isSelected = behavData.task == 1 & behavData.locV ~= behavData.locA;
    behavData.cmb = NaN(size(behavData,1),1);
    behavData.cmb(isSelected) = (behavData.resp(isSelected) - behavData.locA(isSelected)) ...
                    ./(behavData.locV(isSelected) - behavData.locA(isSelected));
    
    % Including to be excluded flag from EEG data
    temp = collecteegdata(subID,expStage,'toReject');
    if ~isempty(temp)
        % Finding the trials in the receiving table which are also
        % present in the table to be appened.
        tempTable = table(temp(:,1),temp(:,2),'VariableNames',{'iSessionOverall','iTrialInSession'});
        [Lia,Locb] = ismember(behavData(:,{'iSessionOverall','iTrialInSession'}),tempTable);
        behavData.toReject = ones(size(behavData,1),1)*9;
        behavData.toReject(Lia) = temp(Locb(Locb ~= 0),3);
    end
    
    % Including prestimulus alpha power from EEG data
    [temp,alphaPowTime] = collecteegdata(subID,expStage,'psalpha');
    if ~isempty(temp)
        % Finding the trials in the receiving table which are also
        % present in the table to be appened.
        tempTable = table(temp(:,1),temp(:,2),'VariableNames',{'iSessionOverall','iTrialInSession'});
        [Lia,Locb] = ismember(behavData(:,{'iSessionOverall','iTrialInSession'}),tempTable);
        behavData.psalpha = NaN(size(behavData,1),size(temp,2)-2);
        behavData.psalpha(Lia,:) = temp(Locb(Locb ~= 0),3:end);
        if isempty(behavData.Properties.UserData)
            behavData.Properties.UserData = ...
                struct('psalpha_cols_time',alphaPowTime);
        else
            behavData.Properties.UserData.('psalpha_cols_time') = ...
                alphaPowTime;
        end
    end
    
    % Saving subject specific behavioural data
    fprintf('\n\nSaving data...\n\n');
    savePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),['preproc_BEHAV_',subID,'.mat']);
    save(savePath,'behavData','-v7.3');
    
    % Organizing data for output
    fieldName = 'behavData';
    subjects(index).subID = str2double(subID);
    subjects(index).(fieldName) = behavData;
    
    index = index + 1;
end

end


function data = collectruns(dirName,type)
% Collects the structs of runs of one subject into a cell array. 
% 
%   Input: 
%       dirName: the full path of the data folder in which a given 
%           subject's data are located. 
%       type: type of files to be used (either 'exp' or 'train')
%
%   Output: cell array of structs, each containing the data of a run. 

currDF = cd(dirName);

listing = dir(['*_',type,'_*_all_present_*.mat']);
fileNames = {listing.name}';

% Sorting the list of file names to natural order. 
fileNames = sort_nat(fileNames); 

% Extracting the structs from the .mat files into a cell array. 
if strcmp(type,'exp')
    
    data = cell(size(fileNames));
    
    for i = 1:numel(fileNames)
        
        load(fullfile(dirName,fileNames{i}),'info');
        
        iDay = repmat(info.iDay,info.inputParameter.nTrials,1);
        iSessionInDay = repmat(info.iSessionInDay,info.inputParameter.nTrials,1);
        iSessionOverall = repmat(info.iSessionOverall,info.inputParameter.nTrials,1);
        iTrialInSession = (1:info.inputParameter.nTrials)';
        condition = info.inputParameter.condition;
        locA = info.inputParameter.locationAuditory;
        locV = info.inputParameter.locationVisual;
        relV = info.inputParameter.reliabilityVisual;
        resp = info.resp.locResp;
        respTime = info.resp.locRespTimes-info.timing.actVisualStartTime;
        task = info.inputParameter.task;
        wrongHand = info.resp.wrongHand;
        
        temp = table(iDay,iSessionInDay,iSessionOverall,iTrialInSession,condition,task,locA,locV,relV,resp,respTime,wrongHand);
                
        data{i} = temp;
        
    end
    
    data = cat(1,data{:});
    data = sortrows(data,{'iDay','iSessionInDay','iTrialInSession'});
    
elseif strcmp(type,'train')
    error('collectruns:unsupportedDataType',...
        'Analysis on training data is not yet implemented!');
end

cd(currDF);

end

function varargout = collecteegdata(subID,expStage,dataToCollect)

switch dataToCollect
    case 'toReject'
        fileMatchStr = 'fteeg_MVPA-sm_.*.mat';
        saveDf = cd(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa_preproc',subID));
        fileList = dir;
        fileList = {fileList.name}';
        matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
    case 'psalpha'
        fileMatchStr = 'fteeg_PSALPHA-tf_[0-9]{3}.mat';
        saveDf = cd(DEC_2_setupdir(expStage,'anal_eeg_sub_psalpha',subID));
        fileList = dir;
        fileList = {fileList.name}';
        matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
end

if sum(matchID) == 0
    warning('collecteegdata:missingData',...
        'No EEG ''%s'' data available for subject %s, skipping! ',...
        dataToCollect,subID);
    cd(saveDf);
    varargout{1} = [];
    return;
end

fileList = fileList(matchID);
data = [];

for i = 1:numel(fileList)
    
    s = load(fileList{i});
    sFieldNames = fieldnames(s);
    s = s.(sFieldNames{1});
    switch dataToCollect
        case 'toReject'
            data = cat(1,data,s.trialinfo(:,[1,2,7]));
        case 'psalpha'
            % Extracting prestim alpha averaged over occipital electrodes
            isAlpha = s.freq >= 8 & s.freq <= 12;
            isOccElec = ismember(s.label,{'PO3','POz','PO4','O1','Oz','O2'});
            temp = s.powspctrm;
            % temp is nTrials x nChannels x nFreq x nTime
            temp = squeeze(nanmean(temp(:,:,isAlpha,:),3));
            % temp is nTrials x nChannels x nTime
            psAlphaOccPowByTrial = squeeze(nanmean(temp(:,isOccElec,:),2));
            % psAlphaOccPowByTrial is nTrials x nTime
            data = cat(1,data,[s.trialinfo(:,[1,2]),psAlphaOccPowByTrial]);
    end
    
end

varargout{1} = data;
varargout{2} = s.time;

cd(saveDf);

end