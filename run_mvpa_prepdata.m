% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

expStage = 'final';

% In case running on server set the path
if ~isempty(regexp(setupID,'^bb.+','once'))
    % Remove unnecessary folders from path preventing collisions
    pathCell = strsplit(path,':')';
    idx = ~cellfun(@isempty,regexp(pathCell,...
                                   '.*ThirdPartyToolboxes.*','once'));
    pathCell = pathCell(idx);
    for i = 1:numel(pathCell), rmpath(pathCell{i}); end
    % Adding necessary folders to path
    addpath(genpath(fullfile(DEC_2_setupdir(expStage,'utils'),'Utility')));
    addpath(fullfile(DEC_2_setupdir(expStage,'utils'),'spm12'));
    addpath(DEC_2_setupdir(expStage,'anal_scripts'));
    spm('defaults', 'EEG');
end

% Input parameters 
% trainMethod = 'sample-wise-source';
trainMethod = 'sample-wise-sm-avg';
subID = {'108','109','110','111','112','113','116','118','119','120','121',...
    '122','123'};
% subID = {'108'};

for i = 1:size(subID,2)
    analysisDir = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID{i}),trainMethod);
    if ~exist(analysisDir,'dir')
        mkdir(analysisDir);
    end
    
    I = struct();
    I.dir_analysis = analysisDir;
    I.dir_preproc = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID{i}),'sample-wise-sm');
    % I.dir_preproc = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa_preproc',subID{i}));
    I.subID = subID{i};
    I.tr_method = trainMethod;
    mvpa.prepdata(I);
end