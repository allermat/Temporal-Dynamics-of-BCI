%% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
    expStage = 'final';
    % Remove unnecessary folders from path preventing collisions
    pathCell = strsplit(path,':')';
    idx = ~cellfun(@isempty,regexp(pathCell,...
                                   '.*ThirdPartyToolboxes.*','once'));
    pathCell = pathCell(idx);
    for i = 1:numel(pathCell), rmpath(pathCell{i}); end
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
        parpool('local',16);
    else
        parpool('local');
    end
end

expStage = 'final';
addpath(genpath(fullfile(DEC_2_setupdir(expStage,'utils'),'Utility')));
fname = 'gr_er_tr-AV-c-av_s_gen-AV-ci-av.mat';

m_ci = mvpares(fullfile(DEC_2_setupdir(expStage,'anal_eeg_group_mvpa'),'sample-wise-sm-avg',fname));
srcFiles = m_ci.getInfo.sourceFiles;

temp = regexp(srcFiles,'(.*)\\(.*\.mat)','tokens','once');
temp = cat(1,temp{:});
p = temp(:,1);
f = temp(:,2);
temp = regexp(p,'\\([0-9]{3})\\','tokens','once');
subIDs = cat(1,temp{:});
fun = @(s,f) fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',s),...
                        'sample-wise-sm-avg',f);
pathIndivFiles = cellfun(fun,subIDs,f,'UniformOutput',false);
% Time window for statistical analysis is based on the time window of
% significant decoding in audiovisual congruent trials. 
m_ci =  mvpa.addStats(m_ci,'avWeights','pathIndivFiles',pathIndivFiles,'timeWin',[0.055,0.7]);
m_ci =  mvpa.addStats(m_ci,'avModelCorr','pathIndivFiles',pathIndivFiles,'timeWin',[0.055,0.7]);
