function saveFigureData()

%% Figure 1
% Behavioural responses
load(fullfile(DEC_2_setupdir('final','anal_behav_group'),'BEHAV_ANAL_group.mat'));
description = ...
    ['Data for Figure 1C, Behavioural Responses. ',...
    'In each row, the structure field ''wAV'' contains the relative audiovisual ',...
    'weight index values observed for each of the thirteen ',...
    'participants (field ''participant'', 1-13), auditory (1) or visual (2) ',...
    'location reports (field ''taskRelevance''), high (1) or low (2) ',...
    'visual reliability (field ''visualReliability''), and high (1) or low (2) ',...
    'audiovisual disparity (field ''disparity'').'];
wAV = grWav.data(grWav.data.groupings == 'all',:);
wAV = wAV{:,3:end};
wAV = wAV(:);
nSubj = numel(unique(grWav.data.subID));
participant = repmat((1:nSubj)',8,1);
visualReliability = repmat([1,1,1,1,2,2,2,2],nSubj,1);
visualReliability = visualReliability(:);
disparity = repmat([1,1,2,2,1,1,2,2],nSubj,1);
disparity = disparity(:);
taskRelevance = repmat([1,2,1,2,1,2,1,2],nSubj,1);
taskRelevance = taskRelevance(:);
fig1Data.behav = struct('description',description,'wAV',wAV,...
                        'participant',participant,...
                        'visualReliability',visualReliability,...
                        'disparity',disparity,...
                        'taskRelevance',taskRelevance);
% Bayesian Causal Inference predictions
load(fullfile(DEC_2_setupdir('final','anal_behav_group'),'bci_simul_wav_BEHAV_group.mat'));
description = ...
    ['Data for Figure 1C, BCI predictions. ',...
    'In each row, the structure field ''wAV'' contains the relative audiovisual ',...
    'weight index values predicted by the BCI model for each of the thirteen ',...
    'participants (field ''participant'', 1-13), auditory (1) or visual (2) ',...
    'location reports (field ''taskRelevance''), high (1) or low (2) ',...
    'visual reliability (field ''visualReliability''), and high (1) or low (2) ',...
    'audiovisual disparity (field ''disparity'').'];
wAV = grWav.data(grWav.data.model == 'bci',:);
wAV = wAV{:,3:end};
wAV = wAV(:);
fig1Data.bciPred = struct('description',description,'wAV',wAV,...
                          'participant',participant,...
                          'visualReliability',visualReliability,...
                          'disparity',disparity,...
                          'taskRelevance',taskRelevance);
save(fullfile(DEC_2_setupdir('final','study_root'),'docs','Data_S1_Figure_1C.mat'),'fig1Data');

%% Figure 2
fig2Data = loadCrossmodGen;
fig2Data.description = ...
    ['The fields ''R_trA_genA'', ''R_trA_genV'', ''R_trV_genA'', ',...
    '''R_trV_genV'' contain the across condition and time generalization ',...
    'matrices for each participant (1-13) trained on unisensory auditory ',...
    '(''trA'') or visual (''trV'') and generalized to unisensory auditory ',...
    '(''genA'') or visual conditions (''genV''). Generalization performance ',...
    'is expressed as the Fisher-transformed correlation coefficient between ',...
    'the true and predicted spatial locations. The dimensionality of these fields is ',...
    'Training time x Generalization time x Participant. Training and ',...
    'generalization time points (in seconds) are given in the fields',...
    '''trTime'' and ''genTime'' respectively'];

save(fullfile(DEC_2_setupdir('final','study_root'),'docs','Data_S1_Figure_2.mat'),'fig2Data');
%% Figure 4
fig4Data.ABCD = loadWav();
fig4Data.ABCD.description = ...
    ['Data for Figure 4A-D ',...
    'Rows of the structure field ''wAV'' contain the relative audiovisual ',...
    'weight index values observed for each of the thirteen ',...
    'participants (field ''participant'', 1-13), auditory (1) or visual (2) ',...
    'location reports (field ''taskRelevance''), high (1) or low (2) ',...
    'visual reliability (field ''visualReliability''), and high (1) or low (2) ',...
    'audiovisual disparity (field ''disparity''). Columns of the field ''wAV'' ',...
    'are the time points (in seconds) as denoted in the field ''time'''];
fig4Data.E = loadWavCorr();
fig4Data.E.description = ...
    ['Data for Figure 4E ',...
    'Rows of the structure field ''corr'' contain the circular-circular ',...
    'correlation coefficients between the behavioural and neural ',...
    'relative audiovisual weight index values observed for each of the thirteen ',...
    'participants (field ''participant'', 1-13). Columns of the field ''corr'' ',...
    'are the time points (in seconds) as denoted in the field ''time'''];
fig4Data.F = loadBCI();
fig4Data.F.description = ...
    ['Data for Figure 4F ',...
    'Rows of the structure field ''bic'' contain the ',...
    'Bayesian Infromation Criterion values from the Bayesien Modelling analysis ',...
    'for each of the five fitted models (field ''models'', 1 - Bayesian ',...
    'Causal Inference, 2 - forced fusion, 3 - segregation unisensory auditory, ',...
    '4 - segregation unisensory visual, 5 - full audiovisual segregation), ',...
    'for each of the thirteen participants (field ''participant'', 1-13). ',...
    'Columns of the field ''bic'' are the time points (in seconds) as ',...
    'denoted in the field ''time'''];

save(fullfile(DEC_2_setupdir('final','study_root'),'docs','Data_S1_Figure_4.mat'),'fig4Data');
end

function data = loadCrossmodGen()

% Collecting subjects
expStage = 'final';
saveDf = cd(DEC_2_setupdir(expStage,'anal_eeg'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);
data.R_trA_genA = [];
data.R_trA_genV = [];
data.R_trV_genA = [];
data.R_trV_genV = [];

for iSubj = 1:numel(subjList)
    subID = subjList{iSubj};
    dir_analysis = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID),...
        'sample-wise-sm-avg');
    matchStr = 'er_tr-[AV]_stim_-100-5-1000_gen-.*';
    
    saveDf = cd(dir_analysis);
    
    fileList = cellstr(ls('*.mat'));
    matchIdx = ~cellfun(@isempty,regexp(fileList,matchStr));
    fileList = fileList(matchIdx);
    
    for i = 1:size(fileList,1)
        switch i
            case 1
                condFname = 'R_trA_genA';
            case 2
                condFname = 'R_trA_genV';
            case 3
                condFname = 'R_trV_genA';
            case 4
                condFname = 'R_trV_genV';
        end
        m = mvpares(fileList{i});
        perf = m.getGenPerfEstimates('genTime','tr_x_tr');
        trTime = m.getTrTimePoints;
        timeMask = trTime >= -0.1 & trTime <= 0.7;
        perf = perf.rf_poolFolds(timeMask,timeMask);
        data.(condFname) = cat(3,data.(condFname),perf);
                
    end
    
    cd(saveDf);
end
[data.trTime,data.genTime] = deal(trTime(timeMask));

end

function data = loadWav()
dataFname = 'gr_er_tr-AV-c-av_s_gen-AV-ci-av.mat';
m = mvpares(fullfile(DEC_2_setupdir('final','anal_eeg_group_mvpa'),...
    'sample-wise-sm-avg',dataFname));
info = m.getInfo;
fileList = info.sourceFiles;
nSubj = numel(fileList);
data.wAV = [];
for i = 1:nSubj
    m = mvpares(fileList{i});
    w = m.getAVmodelWeights('unit','deg','smooth',true);
    fn = fieldnames(w);
    idx = ~cellfun(@isempty,regexp(fn,'\w_\w_\w_wav'));
    fn = sort(fn(idx));
    trTime = m.getTrTimePoints;
    timeMask = trTime >= -0.1 & trTime <= 0.7;
    temp = [];
    for j = 1:numel(fn)
        temp = cat(1,temp,w.(fn{j})(timeMask)');
    end
    data.wAV = cat(1,data.wAV,temp);
end

data.time = trTime(timeMask);
data.participant = repmat(1:nSubj,8,1);
data.participant = data.participant(:);
data.visualReliability = repmat([1,1,1,1,2,2,2,2]',nSubj,1);
data.disparity = repmat([1,1,2,2,1,1,2,2]',nSubj,1);
data.taskRelevance = repmat([1,2,1,2,1,2,1,2]',nSubj,1);

end

function data = loadWavCorr()
dataFname = 'gr_er_tr-AV-c-av_s_gen-AV-ci-av.mat';
m = mvpares(fullfile(DEC_2_setupdir('final','anal_eeg_group_mvpa'),...
    'sample-wise-sm-avg',dataFname));
info = m.getInfo;
fileList = info.sourceFiles;
nSubj = numel(fileList);
data.corr = [];
for i = 1:nSubj
    m = mvpares(fileList{i});
    c = m.getAVmodelCorrelations('behav','smooth',true);
    trTime = m.getTrTimePoints;
    timeMask = trTime >= -0.1 & trTime <= 0.7;
    data.corr = cat(1,data.corr,c.cc_behav(timeMask)');
end

data.time = trTime(timeMask);
data.participant = [1:nSubj]';

end

function data = loadBCI()

% Collecting subjects
expStage = 'final';
saveDf = cd(DEC_2_setupdir(expStage,'anal_eeg'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);
nSubj = numel(subjList);
data.bic = [];
for iSubj = 1:nSubj
    subID = subjList{iSubj};
    dir_analysis = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID),...
        'sample-wise-sm-avg');
    matchStr = ['bci_simul_MVPA_',subID,'_AV-c-av_[0-9]+.mat'];
    
    saveDf = cd(dir_analysis);
    
    fileList = cellstr(ls('*.mat'));
    matchIdx = ~cellfun(@isempty,regexp(fileList,matchStr));
    dataFname = fileList{matchIdx};
    load(dataFname,'bic','timePoints');
    data.bic = cat(1,data.bic,cell2mat(bic)');
    cd(saveDf);
end
data.participant = repmat(1:nSubj,5,1);
data.participant = data.participant(:);
data.model = repmat([1:5]',nSubj,1);
data.time = timePoints';
end