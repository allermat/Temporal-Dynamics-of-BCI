function preproc_EYE(subID)
% DEC_2 project EEG preprocessing common stage

%% Parsing input, checking matlab
p = inputParser;

addRequired(p,'subID',@ischar);

parse(p,subID);

subID = p.Results.subID;

if isempty(regexp(path,'parfor_progress','once'))
    error('parfor_progress not found in path!');
end

%% Opening parallel pool. 
currPool = gcp('nocreate');
% if there is no parallel pool opened, open one.
if isempty(currPool)
    parpool('local');
end

%% Preparing file and directory names for the processing pipeline.
expStage = 'final';
tag = 'EYE';

dataDir = DEC_2_setupdir(expStage,'data_behav_sub',subID);
analysisDir = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub',subID),tag);
if ~exist(analysisDir,'dir')
    mkdir(analysisDir);
end

% Getting data file names
saveDf = cd(dataDir);
fileNames = sort_nat(cellstr(ls([subID,'_exp*eyelink*.asc'])));
cd(saveDf);

% Loading files specifying parameters
trigDef = load('trigDef.mat');
trigDef = trigDef.trigDef;
setupSpec = load('setup_spec.mat');
setupSpec = setupSpec.setup_spec;

stimTriggers = trigDef.trig_eye(trigDef.type == 'stim');
trig_visonset_corr_eyelink = setupSpec.trig_visonset_corr_eyelink;

dataEye = struct('session',[],'event',[],'Fs',[]);
dataEye(size(fileNames,1)).session = [];

fprintf('\n\nConverting files...\n');
parfor_progress(size(fileNames,1));

parfor iFile = 1:size(fileNames,1)
    
    %% Reading raw data
    [event,hdr] = read_eyelink_event(fullfile(dataDir,fileNames{iFile}));
    
    evValues = {event.value}';
    evTypes = {event.type}';
    evStartSamples = [event.sample]';
    % Finding the onset samples of stimulus triggers
    targEvStartSamples = evStartSamples(ismember(evTypes,'Stimulus') & ...
        ismember(evValues,stimTriggers));
    % Correcting stimulus triggers for trigger-visual onset asynchrony
    targEvStartSamples = targEvStartSamples+(trig_visonset_corr_eyelink*hdr.Fs);
    % Saving corrected values into the original structure
    evStartSamples(ismember(evTypes,'Stimulus') & ...
        ismember(evValues,stimTriggers)) = targEvStartSamples;
    evStartSamples = num2cell(evStartSamples);
    [event.sample] = evStartSamples{:};
    
    dataEye(iFile).session = iFile; %#ok<PFOUS>
    dataEye(iFile).event = event;
    dataEye(iFile).Fs = hdr.Fs;
    
    % Advancing Progress monitor
    parfor_progress;
    
end
% Finalizing progress monitor.
parfor_progress(0);

%% Saving data
fprintf('\n\nSaving data...\n\n');
savePath = fullfile(analysisDir,[tag,'_',subID,'.mat']);
save(savePath,'dataEye','-v7.3');
    


%% Closing parallel pool. 
% currPool = gcp('nocreate');
% % If there is a parallel pool opened, close it. 
% if ~isempty(currPool)
%     delete(currPool);
% end

end
