clear all;
clc;
% Collecting subjects
expStage = 'final';
saveDf = cd(DEC_2_setupdir(expStage,'anal_behav'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);


for i = 1:size(subjList,1)
    
    subID = subjList{i};
    
    behavFilePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),...
        ['preproc_BEHAV_',subID,'.mat']);
    if exist(behavFilePath,'file')
        behavData = load(behavFilePath);
        behavData = behavData.dataBehav;
    else
        warning('No file, skipping this subject! ');
        continue;
    end
    
    [behavData.respV,behavData.respA] = deal(NaN(size(behavData,1),1));
    behavData.respV(behavData.task == 2) = behavData.resp(behavData.task == 2);
    behavData.respA(behavData.task == 1) = behavData.resp(behavData.task == 1);
    
    % Selecting appropriate trials and variables for BCI fitting
    bciData = behavData(...
        ~isnan(behavData.resp) & ...    % non-NaN response trials
        ~isnan(behavData.locV) & ...    % Audio-visual trials
        ~isnan(behavData.locA) & ...
        behavData.toReject == 0,...     % Valid trials
        {'locV','locA','relV','respV','respA'});  % and these variables
    
    % Saving subject specific behavioural data
    fprintf('\n\nSaving data...\n\n');
    savePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),...
        ['bci_data_BEHAV_',subID,'.mat']);
    save(savePath,'bciData','-v7.3');
    
end