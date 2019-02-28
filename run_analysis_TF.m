%% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};
expStage = 'final';
if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
%     % Remove unnecessary folders from path preventing collisions
%     pathCell = strsplit(path,':')';
%     idx = ~cellfun(@isempty,regexp(pathCell,...
%                                    '.*ThirdPartyToolboxes.*','once'));
%     pathCell = pathCell(idx);
%     for i = 1:numel(pathCell), rmpath(pathCell{i}); end
    % Adding necessary folders to path
    addpath(genpath(fullfile(DEC_2_setupdir(expStage,'utils'),'Utility')));
    addpath(fullfile(DEC_2_setupdir(expStage,'utils'),'fieldtrip'));
    ft_defaults;
else
    isServer = false;
    dbstop if error;
end

%% Opening parallel pool. 
% if there is no parallel pool running, open one. 
currPool = gcp('nocreate');
if isempty(currPool)
    if isServer
        parpool('local',1);
    else
        parpool('local');
    end
end

% subID = {'108'};
% subID = {'108','109','110','111','112','113','116','118','119','120','121',...
%     '122','123'};
% parfor i = 1:size(subID,2)
%     scratch_TF_baselinecorr(subID{i});
% %     analysis_TF(subID{i});
% end
analysis_TF_group
