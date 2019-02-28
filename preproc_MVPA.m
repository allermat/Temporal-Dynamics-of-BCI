function preproc_MVPA(subID)
%% Parsing input, checking matlab
p = inputParser;

addRequired(p,'subID',@ischar);

parse(p,subID);

subID = p.Results.subID;

if isempty(regexp(path,'parfor_progress','once'))
    error('parfor_progress not found in path!');
end

%% Opening parallel pool. 
% currPool = gcp('nocreate');
% % if there is no parallel pool opened, open one. 
% if isempty(currPool)
%     parpool('local');
% end

%% Preparing file and directory names for the processing pipeline.
expStage = 'final';

eegDataDir = DEC_2_setupdir(expStage,'data_eeg_sub',subID);
sourceDir = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub',subID),'COMM');
behavSourceDir = DEC_2_setupdir(expStage,'anal_behav_sub',subID);
eyeSourceDir = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub',subID),'EYE');
destDir = DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa_preproc',subID);

% Loading files specifying parameters
trigDef = load('trigDef.mat');
trigDef = trigDef.trigDef;
setupSpec = load('setup_spec.mat');
setupSpec = setupSpec.setup_spec;
subjectSpec = load('subject_spec.mat');
subjectSpec = subjectSpec.subject_spec;
% Finding subject specific info
idx = ismember(cellfun(@num2str,{subjectSpec.subID},'UniformOutput',false),subID);
if all(~idx)
    error('No specific info found for this subject.');
end
subjectSpec = subjectSpec(idx);

% Getting original eeg data file names' list
dataFileNames = subjectSpec.eeg_files.fileName;
days = unique(regexp(dataFileNames,'day\d','match','once'));
dataFileNamesInDays = groupListByKey(dataFileNames,'day\d');

% Getting eeg source file names' list
sourceFileNames = strcat('fteeg_COMM_',dataFileNames,'.mat');
% Checking if the required files are available
saveDf = cd(sourceDir);
sourceFileNamesPresent = cellstr(ls(['fteeg_COMM_',subID,'_day*.mat']));
if any(~ismember(sourceFileNames,sourceFileNamesPresent))
    error('The specified source files are not found');
end
sourceFileNamesInDays = groupListByKey(sourceFileNames,'day\d');

% Getting artefact file names' list
artefactFileNames = strcat('artf_COMM_',dataFileNames,'.mat');
% Checking if the required files are available
artefactFileNamesPresent = cellstr(ls(['artf_COMM_',subID,'_day*.mat']));
sourceFileNamesPresent = cellstr(ls(['fteeg_COMM_',subID,'_day*.mat']));
if any(~ismember(artefactFileNames,artefactFileNamesPresent))
    error('The specified artefact files are not found');
end
artefactFileNamesInDays = groupListByKey(artefactFileNames,'day\d');
cd(saveDf);

% Loading behavioural and eyetracking data (these are organized separateli 
% in one single file per subject)
dataEye = load(fullfile(eyeSourceDir,['EYE_',subID,'.mat']));
dataEye = dataEye.dataEye;
dataBehav = load(fullfile(behavSourceDir,['preproc_BEHAV_',subID,'.mat']));
dataBehav = dataBehav.behavData;

% % Keeping track of the already processed sessions, indexing starts from 0
% latestSession = 0;

fprintf('\n\nProcessing files...\n');
parfor_progress(size(sourceFileNames,1));

for iDay = 1:size(dataFileNamesInDays,1)
    
    dataFileNamesActDay = dataFileNamesInDays{iDay};
    sourceFileNamesActDay = sourceFileNamesInDays{iDay};
    artefactFileNamesActDay = artefactFileNamesInDays{iDay};
    
    % Cell array for collecting result files
    resultFilesActDay = cell(numel(sourceFileNamesActDay),1);
    
    for iFileActDay = 1:size(sourceFileNamesActDay,1)
        %% Loading the source data
        ftDataReref = load(fullfile(sourceDir,[sourceFileNamesActDay{iFileActDay},'.mat']));
        ftDataReref = ftDataReref.ftDataReref;
        
        %% low-pass filtering
        cfg = struct();
        cfg.lpfilter = 'yes';
        cfg.lpfreq = 45;
        cfg.lpfiltord = 5;
        cfg.lpfiltdir = 'twopass';
        ftDataLp = ft_preprocessing(cfg,ftDataReref);
        % Clearing unnecessary previous dataset
        ftDataReref = []; 
        
        %% Epoching
        % Trial definition
        cfg = struct();
        cfg.headerfile = fullfile(eegDataDir,[dataFileNamesActDay{iFileActDay},'.vhdr']);
        cfg.trialfun = 'ft_trialfun_eventlocked';
        cfg.trialdef = struct();
        cfg.trialdef.prestim = 0.1;
        cfg.trialdef.poststim = 1;
        cfg.trialdef.trigdef = trigDef;
        cfg.trialdef.fileSpec = subjectSpec.eeg_files(ismember(subjectSpec.eeg_files.fileName,dataFileNamesActDay{iFileActDay}),:);
        cfg.trialdef.trig_visonset_corr = setupSpec.trig_visonset_corr_eeg;
        cfg.trialdef.eventtype = 'Stimulus';
        cfg = ft_definetrial(cfg);
        
        % epoching with background correction
        ftDataEp = ft_redefinetrial(cfg,ftDataLp);
        
        % Clearing unnecessary previous dataset
        ftDataLp = []; 
        
        %% Baseline correction
        cfg = struct();
        cfg.demean = 'yes';
        cfg.baselinewindow = [-0.1,0];
        ftDataEp = ft_preprocessing(cfg,ftDataEp);
        
        %% Rejecting trials based on artefacts, eye events and behaviour. 
                
        % Loading artefact data
        dataArtf = load(fullfile(sourceDir,[artefactFileNamesActDay{iFileActDay},'.mat']));
        
        % Defining behavioural criteria
        % minimum response time in s
        minRespTime = 1;
        
        % Defining saccade thresholds
        % minimum saccade amplitude in degrees
        saccade = struct();
        saccade.ampl = 1;
        % minimum saccade peak velocity in deg/s
        saccade.pv = 15;
        % minimum saccade duration in s
        saccade.dur = 0.06;
        saccadeThresholds = [saccade.ampl,saccade.pv,saccade.dur];
        
        % Defining fixation criteria
        fixation = struct();
        % x coordinate of fixation cross in pixels
        fixation.x = (setupSpec.screen_coords_pix(3)-setupSpec.screen_coords_pix(1))/2;
        % y coordinate of fixation cross in pixels
        fixation.y = (setupSpec.screen_coords_pix(4)-setupSpec.screen_coords_pix(2))/2;
        % radius of the allowed area around the fixation cross in pixels. 
        % (converted from 2 degrees)
        fixation.r = round(posdeg2pix(2,setupSpec.view_dist,setupSpec.screen_coords_pix(3)/setupSpec.screen_coords_mm(3)));
        fixationCriteria = [fixation.x,fixation.y,fixation.r];
        
        % Window of interest around the onset of the stimulus which has to
        % be free of artefacts and unwanted eye events. 
        winOfInt = [-0.1,1];
        
        % Marking bad trials
        % Last column of trialinfo gives the following info:
        % 0 - good
        % 1 - no response
        % 2 - early response
        % 3 - EEG artefact
        % 4 - missing eyetracker data
        % 5 - eyeblink
        % 6 - saccade
        % 7 - wrong fixation location
        badTrials = rejecttrials(ftDataEp,dataBehav,dataEye,dataArtf,trigDef,winOfInt,minRespTime,saccadeThresholds,fixationCriteria,'normal');
        
        ftDataEp.trialinfo = [ftDataEp.trialinfo,badTrials];
        
        resultFilesActDay{iFileActDay} = ftDataEp;
        
        ftDataEp = []; 
        
        % Advancing Progress monitor
        parfor_progress;
        
    end
        
    %% Merging files for the same day if necessary
    if numel(resultFilesActDay) > 1
        ftDataClean = ft_appenddata([],resultFilesActDay{:});
    else
        ftDataClean = resultFilesActDay{1}; %#ok<*NASGU>
    end
    
    resultFilesActDay = [];
    
    %% Saving data
    savePath = fullfile(destDir,['fteeg_MVPA_',subID,'_',days{iDay},'.mat']);
    save(savePath,'ftDataClean','-v7.3');
    
    ftDataClean = [];
    
end

% Finalizing progress monitor.
parfor_progress(0);

% % Closing parallel pool. 
% currPool = gcp('nocreate');
% % If there is a parallel pool opened, close it. 
% if ~isempty(currPool)
%     delete(currPool);
% end

end

function fileNamesByKey = groupListByKey(list,key)

keyLevels = unique(regexp(list,key,'match','once'));

if isempty(keyLevels)
    error('The key is not found in the provided list');
end

fileNamesByKey = cell(numel(keyLevels),1);
for i = 1:numel(keyLevels)
    fileNamesByKey{i} = regexp(list,['\w*',keyLevels{i},'\w*'],'match','once');
    fileNamesByKey{i} = fileNamesByKey{i}(~strcmp(fileNamesByKey{i},''));
end

end