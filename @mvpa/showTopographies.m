function showTopographies(filePath)
% Method for plotting EEG topographies from the MVPA dataset
% 
% DETAILS: further input is given by means of input dialogs during
%   execution
% USAGE:
%   showTopographies(filePath)
% INPUT:
%   filePath (string): path to the MVPA dataset on disk
% OUTPUT: -
%   plots EEG topographies according to the settings given to the input
%       dialogs

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

mvpaData = matfile(filePath);
% info = mvpaData.info;
fieldNames = fieldnames(mvpaData);
featIdx = ~cellfun(@isempty,regexp(fieldNames,'feat'));
featNames = sort(fieldNames(featIdx));
infoIdx = ~cellfun(@isempty,regexp(fieldNames,'info'));
infoNames = sort(fieldNames(infoIdx));
% Select features if there are more than one
if sum(featIdx) > 1
    [s,v] = listdlg('PromptString','Select a feature:','SelectionMode','single',...
                'ListString',featNames);
    if ~v, error('Feature selection failed! '); end
    featToUse = featNames{s};
    info = mvpaData.(infoNames{s});
else
    featToUse = 'feat';
    info = mvpaData.info;
end
% Select session
sessionList = unique(info.session);
sessionListCellStr = cellfun(@num2str,num2cell(sessionList),'UniformOutput',false);
[s,v] = listdlg('PromptString','Select a session:','SelectionMode','single',...
    'ListString',sessionListCellStr);
if ~v, error('Session selection failed! '); end
sessionToUse = sessionList(s);
% Select conditions/response labels
if ~isempty(regexp(featToUse,'_resp','once'))
    respList = info.resp(info.session == sessionToUse);
    respListCellStr = cellfun(@num2str,num2cell(respList),'UniformOutput',false);
    [s,v] = listdlg('PromptString','Select response label:','SelectionMode','multiple',...
        'ListString',respListCellStr);
    if ~v, error('Condition selection failed! '); end
    isExample = info.session == sessionToUse & ismember(info.resp,respList(s));
else
    conditionList = info.condition(info.session == sessionToUse);
    conditionListCellStr = cellfun(@num2str,num2cell(conditionList),'UniformOutput',false);
    [s,v] = listdlg('PromptString','Select conditions:','SelectionMode','multiple',...
        'ListString',conditionListCellStr);
    if ~v, error('Condition selection failed! '); end
    isExample = info.session == sessionToUse & ismember(info.condition,conditionList(s));
end

misc = mvpaData.misc;
[~,nTimePoints,~] = size(mvpaData,featToUse);
time = (0:nTimePoints-1/misc.fs)+misc.timeOnset;
if isfield(misc,'feature_labels');
    feature_labels = misc.feature_labels;
else
    feature_labels = {'Fp1';'Fp2';'F7';'F3';'Fz';'F4';'F8';'FC5';'FC1';...
        'FC2';'FC6';'T7';'C3';'Cz';'C4';'T8';'TP9';'CP5';'CP1';'CP2';...
        'CP6';'TP10';'P7';'P3';'Pz';'P4';'P8';'PO9';'O1';'Oz';'O2';...
        'PO10';'AF7';'AF3';'AF4';'AF8';'F5';'F1';'F2';'F6';'FT9';'FT7';...
        'FC3';'FC4';'FT8';'FT10';'C5';'C1';'C2';'C6';'TP7';'CP3';'CPz';...
        'CP4';'TP8';'P5';'P1';'P2';'P6';'PO7';'PO3';'POz';'PO4';'PO8'};
end
exampleIdx = find(isExample);
ftDataStructs = cell(1,sum(isExample));
% EEG data are converted to fieldtrip's 'timelock' data structure. See the
% documentation of ft_datatype_timelock for further information. 
for i = 1:sum(isExample)
    ftDataStructs{i}.dimord = 'chan_time';
    ftDataStructs{i}.avg = mvpaData.(featToUse)(:,:,exampleIdx(i));
    ftDataStructs{i}.label = feature_labels;
    ftDataStructs{i}.time = time;
end

cfg = struct();
cfg.layout       = 'EEG1010.lay';
cfg.interactive  = 'yes';
cfg.colorbar     = 'yes';
cfg.marker       = 'labels';
figure(), ft_multiplotER(cfg,ftDataStructs{:});

end