function preproc_electrode_pos(subID)
%% Parsing input, checking matlab
p = inputParser;

addRequired(p,'subID',@ischar);

parse(p,subID);

subID = p.Results.subID;

% Preparing file and directory names for the processing pipeline.
expStage = 'final';

eegDataDir = DEC_2_setupdir(expStage,'data_eeg_sub',subID);
destDir = DEC_2_setupdir(expStage,'anal_eeg_sub_sourceloc',subID);

% Loading files specifying parameters
subjectSpec = load('subject_spec.mat');
subjectSpec = subjectSpec.subject_spec;
% Finding subject specific info
idx = ismember(cellfun(@num2str,{subjectSpec.subID},'UniformOutput',false),subID);
if all(~idx)
    error('No specific info found for this subject.');
end
subjectSpec = subjectSpec(idx);

% loadning channel map
channelMap = load('actiCAP_64Ch_channel_map.mat');
channelMap = channelMap.actiCAP_64Ch_channel_map;

% Loading the structural mri
mri = load(fullfile(destDir,'mri_realigned'));
mri = mri.mri_realigned;

% Getting the positions of the fiducials
nas = mri.cfg.fiducial.nas;
lpa = mri.cfg.fiducial.lpa;
rpa = mri.cfg.fiducial.rpa;
 
transm = mri.transform;
 
nas = ft_warp_apply(transm,nas, 'homogenous');
lpa = ft_warp_apply(transm,lpa, 'homogenous');
rpa = ft_warp_apply(transm,rpa, 'homogenous');

% create a structure similar to a template set of electrodes
fid.elecpos = [nas; lpa; rpa];
% Fiducial names in elec
fid.label = {'Nasion','LPA','RPA'};
fid.unit = 'mm';

% Getting eeg channel position file names' list
electrodePosFileNames = subjectSpec.eeg_electrode_pos_files.fileName;
days = subjectSpec.eeg_electrode_pos_files.day;

for i = 1:numel(electrodePosFileNames)

    % Loading chanel position file
    elec = ft_read_sens(fullfile(eegDataDir,electrodePosFileNames{i}));
    elec = ft_convert_units(elec, 'mm');
    % Replacing the number labels in the electrode definition with
    % channel names
    elec.label(1:64) = channelMap(1:64,3);
    % Making sure the channel positions are identical to the
    % electrode positions, otherwise the alignment step does not work
    elec.chanpos = elec.elecpos;
    
    % Align the electrodes to the fiducials from the mri
    cfg = [];
    cfg.method = 'fiducial';
    cfg.target = fid;
    cfg.elec = elec;
    cfg.fiducial = {'Nasion', 'LPA', 'RPA'};
    elec_aligned = ft_electroderealign(cfg); %#ok
    
    savePath = fullfile(destDir,['electrode_pos_',subID,'_day',num2str(days(i)),'.mat']);
    save(savePath,'elec_aligned');
    
end


end