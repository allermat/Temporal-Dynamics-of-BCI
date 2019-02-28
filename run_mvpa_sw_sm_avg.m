clc;
clear all;

%% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
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

%% Adding the necessary folders to the path
expStage = 'final';

addpath(fullfile(DEC_2_setupdir(expStage,'utils'),'libsvm-3.20','matlab'));
addpath(genpath(fullfile(DEC_2_setupdir(expStage,'utils'),'Utility')));

%% Input parameters 
subID = {'108','109','110','111','112','113','116','118','119','120','121',...
    '122','123'};

trainMethod = 'sample-wise-sm-avg';

% Settings for training
trSettings(1).cond = 'A';
trSettings(1).label = 'stim';
trSettings(1).fileMatchStr = '';
trSettings(2).cond = 'V';
trSettings(2).label = 'stim';
trSettings(2).fileName = '';
trSettings(3).cond = 'AV-c-av';
trSettings(3).label = 'stim';
trSettings(3).fileName = '';
% Settings for generalizing
genSettings(1).cond = {'A','V'};
genSettings(1).genTime = {'tr_x_tr','tr_x_tr'};
genSettings(2).cond = {'A','V'};
genSettings(2).genTime = {'tr_x_tr','tr_x_tr'};
genSettings(3).cond = {'AV-c-av','AV-ci-av'};
genSettings(3).genTime = {'tr','tr'};

% Whether to display progress monitor
if isServer, progrMonitor = false; else progrMonitor = true; end

for iSubj = 1:numel(subID)
    
    dataFileName = [subID{iSubj},'_sw-sm-avg_data.mat'];
    
    for iTrain = 1:size(trSettings,2)
        
        if isempty(trSettings(iTrain).fileMatchStr)
            I = struct();
            I.cv_scheme = 'loso';
            I.dir_analysis = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID{iSubj}),trainMethod);
            I.dir_dataFile = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID{iSubj}),trainMethod,dataFileName);
            I.progrMonitor = progrMonitor;
            I.sc_method = 'Z-e';
            I.subID = subID{iSubj};
            I.svm_type = 'er';
            I.tr_cond = trSettings(iTrain).cond;
            I.tr_label = trSettings(iTrain).label;
            I.tr_method = trainMethod;
            I.tr_timePoints = num2cell(-100:5:1000);
            
            [~,fn] = mvpa.svmTrain(I);
        end
        
        if ~isempty(genSettings(iTrain))
            
            for iGen = 1:numel(genSettings(iTrain).cond)
                I = struct();
                if isempty(trSettings(iTrain).fileMatchStr)
                    I.dir_trainFile = fn;
                else
                    saveDf = cd(fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID{iSubj}),...
                        trainMethod));
                    listing = dir;
                    fileNames = {listing.name}';
                    fileNames = regexp(fileNames,trSettings(iTrain).fileMatchStr,'match','once');
                    fileNames = fileNames(~strcmp(fileNames,''));
                    if size(fileNames,1) > 1
                        error('More than one file found with the specified match string!')
                    else
                        fileName = fileNames{1};
                    end
                    cd(saveDf);
                    I.dir_trainFile = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID{iSubj}),...
                        trainMethod,fileName);
                end
                I.dir_dataFile = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID{iSubj}),...
                    trainMethod,dataFileName);
                I.gen_cond = genSettings(iTrain).cond{iGen};
                I.gen_time = genSettings(iTrain).genTime{iGen};
                I.progrMonitor = progrMonitor;
                
                mvparesObj = mvpa.svmGeneralize(I);
                if strcmp(I.gen_cond,'AV-ci-av')
                    mvparesObj = mvpa.addAVmodelEstimates(mvparesObj);
                else
                    mvparesObj = mvpa.addGenPerfEstimates(mvparesObj);
                end
            end
            
        end
        
    end
    
end

if ~isServer, dbclear if error; end