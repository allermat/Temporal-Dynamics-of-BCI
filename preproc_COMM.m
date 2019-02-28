function preproc_COMM(subID)
% DEC_2 project EEG preprocessing common stage

%% Parsing input, checking matlab
p = inputParser;

addRequired(p,'subID',@ischar);

parse(p,subID);

subID = p.Results.subID;

%% Preparing file and directory names for the processing pipeline.
expStage = 'final';

dataDir = DEC_2_setupdir(expStage,'data_eeg_sub',subID);
if ~exist(dataDir,'dir')
    error('The specified data folder does not exist!');
end

analysisDir = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub',subID),'COMM');
if ~exist(analysisDir,'dir');
    mkdir(analysisDir);
end

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

% Getting data file names
fileNames = subjectSpec.eeg_files.fileName;
days = unique(regexp(fileNames,'day\d','match','once'));

sourceFileNamesInDays = cell(numel(days),1);
for iDay = 1:numel(days)
    sourceFileNamesInDays{iDay} = regexp(fileNames,['\w*',days{iDay},'\w*'],'match','once');
    sourceFileNamesInDays{iDay} = sourceFileNamesInDays{iDay}(~strcmp(sourceFileNamesInDays{iDay},''));
end

tag = 'COMM_';

for iDay = 1:numel(days)
    
    sourceFileNamesActDay = sourceFileNamesInDays{iDay};
    
    for iFileActDay = 1:size(sourceFileNamesActDay,1)
        
        %% Checking if the processed data are already present
        artfResultFileName = ['artf_',tag,sourceFileNamesActDay{iFileActDay},'.mat'];
        eegResultFileName = ['fteeg_',tag,sourceFileNamesActDay{iFileActDay},'.mat'];
        if exist(fullfile(analysisDir,artfResultFileName),'file') && exist(fullfile(analysisDir,eegResultFileName),'file')
            warning('Skipping %s as it has been already processed! ',sourceFileNamesActDay{iFileActDay});
            continue;
        end
        
        %% Reading raw data
        cfg = struct();
        cfg.datafile = fullfile(dataDir,[sourceFileNamesActDay{iFileActDay},'.eeg']);
        cfg.headerfile = fullfile(dataDir,[sourceFileNamesActDay{iFileActDay},'.vhdr']);
        cfg.channel = 'all';
        ftDataRaw = ft_preprocessing(cfg);
        
        %% High-pass filtering
        cfg = struct();
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 0.1;
        cfg.hpfiltord = 4;
        cfg.hpfiltdir = 'twopass';
        cfg.plotfiltresp = 'yes';
        fprintf('\nFILTERING: Highpass filter of order %d, half power frequency %d\n\n',cfg.hpfiltord,cfg.hpfreq);
        ftDataHp = ft_preprocessing(cfg,ftDataRaw);
        
        %% Re-referencing to average reference
        cfg = struct();
        cfg.reref = 'yes';
        cfg.refchannel = 'all';
        ftDataReref = ft_preprocessing(cfg,ftDataHp);
        % Clearing unnecessary previous dataset
        ftDataHp = []; %#ok<NASGU>
        
        %% Artefact detection and visual check
        if ~exist(fullfile(analysisDir,artfResultFileName),'file')
            %% Automatic artefact detection
            % Marking trials for artefact detection
            cfg = struct();
            cfg.headerfile = fullfile(dataDir,[sourceFileNamesActDay{iFileActDay},'.vhdr']);
            cfg.trialdef.eventtype = 'Stimulus';
            cfg.trialdef.faketrllength = 15;
            cfg.trialdef.prestim = 0.1;
            cfg.trialdef.poststim = 1;
            cfg.trialdef.trigdef = trigDef;
            cfg.trialdef.fileSpec = subjectSpec.eeg_files(ismember(subjectSpec.eeg_files.fileName,sourceFileNamesActDay{iFileActDay}),:);
            cfg.trialdef.trig_visonset_corr = setupSpec.trig_visonset_corr_eeg;
            cfg.trialfun = 'ft_trialfun_artefactdetection';
            cfg = ft_definetrial(cfg);
            trlArtf = cfg.trl;
            
            % Epoching data for automatic artefact detection
            cfg = struct();
            cfg.trl = trlArtf;
            ftDataEpArtf = ft_redefinetrial(cfg,ftDataReref);
            
            % - MUSCLE ARTEFACTS -
            cfg = struct();
            cfg.continuous = 'no';
            % channel selection, cutoff and padding
            cfg.artfctdef.zvalue.channel = 'all'; %{'FT7','FT8','AF7','AF8','TP7','TP8','PO7','PO8'};
            cfg.artfctdef.zvalue.cutoff = subjectSpec.cutoff_zval;
            cfg.artfctdef.zvalue.trlpadding = 0;
            cfg.artfctdef.zvalue.fltpadding = 0;
            cfg.artfctdef.zvalue.artpadding = 0.1;
            % algorithmic parameters
            cfg.artfctdef.zvalue.bpfilter = 'yes';
            cfg.artfctdef.zvalue.bpfreq = [110 140];
            cfg.artfctdef.zvalue.bpfilttype = 'but';
            cfg.artfctdef.zvalue.bpfiltord = 9;
            cfg.artfctdef.zvalue.hilbert = 'yes';
            cfg.artfctdef.zvalue.boxcar = 0.5;
            % make the process interactive
            cfg.artfctdef.zvalue.interactive = 'no';
            
            [~,artefact_muscle] = ft_artifact_zvalue(cfg,ftDataEpArtf);
            
            % - JUMP ARTEFACTS -
            %         cfg = struct();
            %         % required fields
            %         cfg.continuous                  = 'no';
            %         cfg.artfctdef.zvalue.channel    = 'all';
            %         cfg.artfctdef.zvalue.cutoff     = 60; % µV (let's assume that if there are jumps, they are pretty big)
            %         cfg.artfctdef.zvalue.trlpadding = 0;  % add a bit of data on both sides
            %         cfg.artfctdef.zvalue.fltpadding = 0;  % for trial and filter padding
            %         cfg.artfctdef.zvalue.artpadding = 0;  % delete time period on both sides
            %         % algorithmic parameters
            %         cfg.artfctdef.zvalue.cumulative    = 'yes';
            %         cfg.artfctdef.zvalue.medianfilter  = 'yes'; % preserves jumps
            %         cfg.artfctdef.zvalue.medianfiltord = 9;
            %         cfg.artfctdef.zvalue.absdiff       = 'yes';
            %         cfg.artfctdef.zvalue.interactive   = 'no';
            %
            %         [~,artefact_jump] = ft_artifact_zvalue(cfg,dataEpArtf);
            
            %% Inspect detected artifacts
            cfg = struct();
            cfg.viewmode  = 'vertical';
            cfg.continuous = 'no';
            cfg.channel = 'all';
            cfg.artfctdef.muscle.artifact = artefact_muscle;
            %         cfg.artfctdef.jump.artifact = artefact_jump;
            cfg = ft_databrowser(cfg,ftDataEpArtf);
            
            %% Saving artefact data
            fprintf('\n\nSaving artefact data...\n\n');
            % artefacts
            artefact_muscle = cfg.artfctdef.muscle.artifact; %#ok<NASGU>
            artefact_visual = cfg.artfctdef.visual.artifact; %#ok<NASGU>
            save(fullfile(analysisDir,artfResultFileName),'artefact_muscle','artefact_visual','-v7.3');
        end
        
        %% Saving EEG data
        fprintf('\n\nSaving EEG data...\n\n');
        save(fullfile(analysisDir,eegResultFileName),'ftDataReref','-v7.3');
        
    end
end

end
