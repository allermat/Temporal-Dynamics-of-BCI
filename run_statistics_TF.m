%% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
    expStage = 'final';
    % Remove unnecessary folders from path preventing collisions
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
        parpool('local',8);
    else
        parpool('local');
    end
end

% expDesign = {'3way','4way'};
% matchStr = {'(r[hl]-d[hl]-[av])','(r[hl]-[av]-v[lr]-a[lr])'};
expDesign = {'3way','oneSample'};
matchStr = {'(r[hl]-d[hl]-[av])','(r[hl]-d[hl]-[av])'};
for i = 2 %1:numel(expDesign)
    % Low-frequency range
    % Load files
    expStage = 'final';
    sourceDf = DEC_2_setupdir(expStage,'anal_eeg_group_tf');
    saveDf = cd(sourceDf);
    listing = dir;
    condStr = regexp({listing.name}',...
        ['fteeg_TF-pow-low_bis-',matchStr{i},'_gravg.mat'],'tokens','once');
    condStr = strrep([condStr{:}],'-','_');
    listing = regexp({listing.name}',...
        ['fteeg_TF-pow-low_bis-',matchStr{i},'_gravg.mat'],'match','once');
    fNames = listing(~strcmp(listing,''));
    ftData = cell2struct(struct2cell(cellfun(@load,fNames)),condStr,2);
    cd(saveDf);
    % perform stats
    stats = compEEGstats(ftData,expDesign{i},'freq');
    % save data
    save(fullfile(sourceDf,['stat_TF-pow-low_',expDesign{i},'.mat']),...
        'stats','-v7.3');
    
    % High-frequency range
    expStage = 'final';
    sourceDf = DEC_2_setupdir(expStage,'anal_eeg_group_tf');
    saveDf = cd(sourceDf);
    listing = dir;
    condStr = regexp({listing.name}',...
        ['fteeg_TF-pow-high_bis-',matchStr{i},'_gravg.mat'],'tokens','once');
    condStr = strrep([condStr{:}],'-','_');
    listing = regexp({listing.name}',...
        ['fteeg_TF-pow-high_bis-',matchStr{i},'_gravg.mat'],'match','once');
    fNames = listing(~strcmp(listing,''));
    ftData = cell2struct(struct2cell(cellfun(@load,fNames)),condStr,2);
    cd(saveDf);
    stats = compEEGstats(ftData,expDesign{i},'freq');
    save(fullfile(sourceDf,['stat_TF-pow-high_',expDesign{i},'.mat']),...
        'stats','-v7.3');
end