function prepdata(I)
% Prepares EEG data for mvpa. 
%
% INPUT: strucure I with the following fields   
%   dir_analysis:
%   dir_condDefFile:
%   dir_preproc:
%   
%   subID: 
%   tr_method: sample-wise, sample-wise-sm,...


%% Parsing input
p = inputParser;

requiredVars = {'dir_analysis','dir_preproc','subID','tr_method'};
validTrMethods = {'sample-wise','sample-wise-avg','sample-wise-beta',...
    'sample-wise-beta-rl','sample-wise-bp','sample-wise-bp-avg',...
    'sample-wise-rl','sample-wise-sm','sample-wise-sm-corr1',...
    'sample-wise-sm-corr2','sample-wise-sm-avg-corr2','sample-wise-sm-avg',...
    'sample-wise-sm-beta','sample-wise-sm-pseudoav','sample-wise-source',...
    'sample-wise-source-avg','sample-wise-tf'};

addParameter(p,'dir_analysis','',@(x)exist(x,'dir'));
addParameter(p,'dir_preproc','',@(x)exist(x,'dir'));
addParameter(p,'subID','',@ischar);
addParameter(p,'tr_method','',@(x)any(validatestring(x,validTrMethods)));

parse(p,I);

if any(ismember(requiredVars,p.UsingDefaults))
    error('mvpa:prepdata:missingInput',...
        'All required parameters must be specified!');
end

dirAnalysis = p.Results.dir_analysis;
dirPreproc = p.Results.dir_preproc;
subID = p.Results.subID;
trMethod = p.Results.tr_method;

%% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
else
    isServer = false;
end

%%
% Load condition definition file
condDef = load('condDef.mat');
condDef = condDef.condDef;

%% Extracting data from the available files
if strcmp(trMethod,'sample-wise')
    
    destFname = [subID,'_','sw','_','data.mat'];
    eegFileIDstr = 'fteeg_MVPA_';
    
    [feat,info,misc] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod); %#ok<*NASGU,*ASGLU>
    
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'condDef','feat','info','misc','-v7.3');

elseif strcmp(trMethod,'sample-wise-sm')
    
    dataFname = [subID,'_','sw-sm','_','data.mat'];
    eegFileIDstr = 'fteeg_MVPA-sm_';
    
    [feat,info,misc] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod);
    
    % Smoothing the time courses of each trial
    % Window size of moving average, since the sampling rate is 1000 Hz,
    % sample numbers directly correspond to ms.
    mavgwin = 20;
    % Kernel for smoothing trials with the given time wintow
    kern = ones(1,mavgwin)./mavgwin;
    
    % if there is no parallel pool running, open one.
    if isServer
        currPool = gcp('nocreate');
        if isempty(currPool)
            parpool('local',16);
        end
    else
        currPool = gcp('nocreate');
        if isempty(currPool)
            parpool('local');
        end
    end
    
    for i = 1:size(feat,3)
        
        feattmp = feat(:,:,i);
        featsmtmp = feat(:,:,i);
        
        parfor j = 1:size(feattmp,1)
            featsmtmp(j,:) = conv(feattmp(j,:),kern,'same');
        end
        
        feat(:,:,i) = featsmtmp;
    end
    
    % Saving the data
    save(fullfile(dirAnalysis,dataFname),'condDef','feat','info','misc','-v7.3');
    
elseif strcmp(trMethod,'sample-wise-rl')
    
    destFname = [subID,'_','sw-rl','_','data.mat'];
    eegFileIDstr = 'fteeg_MVPA-rl_';
    
    [feat,info,misc] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod);
    
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'condDef','feat','info','misc','-v7.3');
    
elseif any(strcmp(trMethod,{'sample-wise-avg','sample-wise-sm-avg','sample-wise-bp-avg'}))
    
    if strcmp(trMethod,'sample-wise-avg')
        sourceFname = [subID,'_','sw','_','data.mat'];
        destFname = [subID,'_','sw-avg','_','data.mat'];
    elseif strcmp(trMethod,'sample-wise-sm-avg')
        sourceFname = [subID,'_','sw-sm','_','data.mat'];
        destFname = [subID,'_','sw-sm-avg','_','data.mat'];
    elseif strcmp(trMethod,'sample-wise-bp-avg')
        sourceFname = [subID,'_','sw-bp','_','data.mat'];
        destFname = [subID,'_','sw-bp-avg','_','data.mat'];
    end
    M = matfile(fullfile(dirPreproc,sourceFname));
    info = M.info;
    feat = M.feat;
    misc = M.misc;
    
    [info,infoAvg,featAvg] = averagetrials(info,feat,condDef);
    % Preparing variables for saving
    misc.info_source = info;
    info = infoAvg;
    feat = featAvg;
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'condDef','feat','info','misc','-v7.3');

elseif any(strcmp(trMethod,{'sample-wise-sm-corr1'}))
    
    sourceFname = [subID,'_','sw-sm','_','data.mat'];
    destFname = [subID,'_','sw-sm-corr1','_','data.mat'];
    
    M = matfile(fullfile(dirPreproc,sourceFname));
    info = M.info;
    feat = M.feat;
    misc = M.misc;
    
    % Method for generating pseudo-AV data
    [info,infoCorr,featCorr] = mvpa.generateCorrectedAvData(info,feat,condDef,1,'single_trial');
    % Preparing variables for saving
    misc.info_source = info;
    info = infoCorr;
    feat = featCorr;
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'condDef','feat','info','misc','-v7.3');

elseif any(strcmp(trMethod,{'sample-wise-sm-corr2'}))
    
    sourceFname = [subID,'_','sw-sm','_','data.mat'];
    destFname = [subID,'_','sw-sm-corr2','_','data.mat'];
    
    M = matfile(fullfile(dirPreproc,sourceFname));
    info = M.info;
    feat = M.feat;
    misc = M.misc;
    
    % Method for generating pseudo-AV data
    [info,infoCorr,featCorr] = mvpa.generateCorrectedAvData(info,feat,condDef,2,'single_trial');
    % Preparing variables for saving
    misc.info_source = info;
    info = infoCorr;
    feat = featCorr;
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'condDef','feat','info','misc','-v7.3');

elseif any(strcmp(trMethod,{'sample-wise-sm-avg-corr2'}))
    
    sourceFname = [subID,'_','sw-sm-avg','_','data.mat'];
    destFname = [subID,'_','sw-sm-avg-corr2','_','data.mat'];
    
    M = matfile(fullfile(dirPreproc,sourceFname));
    info = M.info;
    feat = M.feat;
    misc = M.misc;
    
    % Method for generating pseudo-AV data
    [info,infoCorr,featCorr] = mvpa.generateCorrectedAvData(info,feat,condDef,2,'avg_potential');
    % Preparing variables for saving
    misc.info_source = info;
    info = infoCorr;
    feat = featCorr;
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'condDef','feat','info','misc','-v7.3');

elseif any(strcmp(trMethod,{'sample-wise-sm-pseudoav'}))
    
    sourceFname = [subID,'_','sw-sm','_','data.mat'];
    destFname = [subID,'_','sw-sm-pseudoav','_','data.mat'];
    
    M = matfile(fullfile(dirPreproc,sourceFname));
    info = M.info;
    feat = M.feat;
    misc = M.misc;
    
    % Method for generating pseudo-AV data
    [info,infoPseu,featPseu] = mvpa.generatePseudoAvData(info,feat,condDef);
    % Preparing variables for saving
    misc.info_source = info;
    info = infoPseu;
    feat = featPseu;
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'condDef','feat','info','misc','-v7.3');

elseif strcmp(trMethod,'sample-wise-source')
    
    destFname = [subID,'_','sw-source','_','data.mat'];
    % V1-3
    eegFileIDstr = 'cspmeeg_MVPA_SOURCE_[0-9]{3}_V1-3.mat';
    [feat_V1_3,info,misc] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod);
    % Saving the first part of the file
    save(fullfile(dirAnalysis,destFname),'condDef','feat_V1_3','info','misc','-v7.3');
    % Cleaning up memory from the big variable(s)
    feat_V1_3 = [];
    % A
    eegFileIDstr = 'cspmeeg_MVPA_SOURCE_[0-9]{3}_A.mat';
    [feat_A,~,~] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod);
    save(fullfile(dirAnalysis,destFname),'feat_A','-append');
    % Cleaning up memory from the big variable(s)
    feat_A = [];
    % IPS
    eegFileIDstr = 'cspmeeg_MVPA_SOURCE_[0-9]{3}_IPS.mat';
    [feat_IPS,~,~] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod);
    save(fullfile(dirAnalysis,destFname),'feat_IPS','-append');
    % feat_* fields will be 3 dimensional: 
    % 1. number of vertices (sources)
    % 2. number of time samples
    % 3. number of examples

elseif ismember(trMethod,{'sample-wise-source-avg'})
    
    sourceFname = [subID,'_sw-source_data.mat'];
    destFname = [subID,'_sw-source-avg_data.mat'];
    
    M = matfile(fullfile(dirPreproc,sourceFname));
    % Getting feature fields for ROIs
    fieldNames = fieldnames(M);
    featIdx = ~cellfun(@isempty,regexp(fieldNames,'feat_.*','once'));
    featFieldNames = fieldNames(featIdx);
    info = M.info;
    
    saveStruct = struct();
    saveStruct.misc = M.misc;
    for i = 1:numel(featFieldNames)
        feat = M.(featFieldNames{i});
        if i == 1
            [info,infoAvg,featAvg] = averagetrials(info,feat,condDef);
        else
            % This way the same trials are averaged consistently
            % across ROIs
            [~,~,featAvg] = averagetrials(info,feat,condDef,'reuseInfo',true);
        end
        
        % Preparing variables for saving
        saveStruct.misc.info_source = info;
        saveStruct.info = infoAvg;
        saveStruct.(featFieldNames{i}) = featAvg;
        saveStruct.condDef = condDef;
    end
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'-struct','saveStruct','-v7.3');
    
elseif any(strcmp(trMethod,{'sample-wise-beta','sample-wise-beta-rl','sample-wise-sm-beta'}))
    
    if strcmp(trMethod,'sample-wise-beta')
        sourceFname = [subID,'_','sw','_','data.mat'];
        destFname = [subID,'_','sw-beta','_','data.mat'];
    elseif strcmp(trMethod,'sample-wise-beta-rl')
        sourceFname = [subID,'_','sw-rl','_','data.mat'];
        destFname = [subID,'_','sw-beta-rl','_','data.mat'];
    elseif strcmp(trMethod,'sample-wise-sm-beta')
        sourceFname = [subID,'_','sw-sm','_','data.mat'];
        destFname = [subID,'_','sw-sm-beta','_','data.mat'];
    end
    M = matfile(fullfile(dirPreproc,sourceFname));
    info = M.info;
    rawData = M.feat;
    misc = M.misc;
    % Finding the number of sessions overall
    nSessionAll = max(info.session);
    % Nubmer of sessions in a day (ugly,ugly hard-wiring...)
    nSessionInDay = 20;
    % Number of days 
    nDays = nSessionAll/nSessionInDay;
    
    % Finding AV condition labels
    AVconds = condDef.condition(~isnan(condDef.locationAuditory) & ~isnan(condDef.locationVisual));
    % Selecting AV trials from the dataset
    AVtrials = ismember(info.condition,AVconds);
    info = info(AVtrials,:);
    rawData = rawData(:,:,AVtrials);
    
    % Defining new sessions by concatenating several sessions
    % Number of new sessions for the dataset
    nNewSession = 12;
    % The concatenation happens separately for each day to reduce the
    % variability of data in a given new session. 
    info.day = NaN(size(info,1),1);
    for iDay = 1:nDays
        startIdx = ((iDay-1)*nSessionInDay)+1;
        endIdx = ((iDay-1)*nSessionInDay)+nSessionInDay;
        info.day(ismember(info.session,startIdx:endIdx)) = iDay;
    end
    % Assigning new sessions to old sessions
    newSessionAssignment = assignnewsessions(info,nNewSession);
    info.newSession = info.session;
    for i = 1:size(newSessionAssignment,1)
        info.newSession(info.session == newSessionAssignment.session(i)) = newSessionAssignment.newSession(i);
    end
    % Defining new variables to make computations easier. 
    info.lowVisRel = info.relV == max(info.relV);
    info.highVisRel = info.relV == min(info.relV);
    info.absAVdiscr = abs(info.locV - info.locA);
    info.smallDiscr = info.absAVdiscr <= 6.67;
    info.largeDiscr = info.absAVdiscr > 6.67;
    % Openin parallel pool to reduce computation time
    if isServer
        currPool = gcp('nocreate');
        if isempty(currPool)
            parpool('local',16);
        end
    else
        currPool = gcp('nocreate');
        if isempty(currPool)
            parpool('local');
        end
    end
    % Whether to include the intercept term
    intercept = false;
    % Computing beta images
    [feat_stim,feat_resp,info_stim,info_resp,design_stim,design_resp] = compbetaimages(info,rawData,intercept,~isServer);
    % Saving additional misc data
    misc.info_source = info;
    misc.design_stim = design_stim;
    misc.design_resp = design_resp;
    % Saving the data
    save(fullfile(dirAnalysis,destFname),'condDef','feat_stim','info_stim','feat_resp','info_resp','misc','-v7.3');
    
elseif strcmp(trMethod,'sample-wise-bp')
    
    destFname = [subID,'_','sw-bp','_','data.mat'];
    
    freqStr = {'d','t','a','b','gl','gh'};
    
    % Loading the first part of the dataset to determine the size. 
    eegFileIDstr = ['fteeg_MVPA-bp-',freqStr{1},'_'];
    [temp,~,~] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod);
    tempSize = size(temp);
    temp = [];
    % Pre-allocating and saving the final array
    feat = NaN(tempSize(1),size(freqStr,2),tempSize(2),tempSize(3));
    save(fullfile(dirAnalysis,destFname),'feat','-v7.3');
    feat = [];
    % Opening the final array in writable mode
    M = matfile(fullfile(dirAnalysis,destFname),'Writable',true);
    % Pre-allocating array for misc information
    misc = cell(size(freqStr));
    for j = 1:numel(freqStr)
        eegFileIDstr = ['fteeg_MVPA-bp-',freqStr{j},'_'];
        [temp,info,misc{j}] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod);
        % Saving the actual part of the array to disk
        M.feat(:,j,:,:) = reshape(temp,size(temp,1),1,size(temp,2),size(temp,3));
        temp = [];
    end
    M = [];
    % feat will be 4 dimensional:
    % 1. number of channels
    % 2. number of frequencies
    % 3. number of time samples
    % 4. number of examples
    temp = [misc{:}];
    misc = temp(1);
    misc.frequencies = {temp.frequencies};
    misc.freqStr = freqStr;
    
    % Appending the remaining variables
    save(fullfile(dirAnalysis,destFname),'condDef','info','misc','-append');
    
elseif strcmp(trMethod,'sample-wise-tf')
    
    destFname = [subID,'_','sw-tf','_','data.mat'];
    % Instantaneous power files
    eegFileIDstr = 'fteeg_MVPA_TF';
    [feat_f,info,misc] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod);
    % Saving the first part of the file
    save(fullfile(dirAnalysis,destFname),'condDef','feat_f','info','misc','-v7.3');
    % Cleanin up memory from the big variable(s)
    feat_f = [];
    % Instantaneous phase files
    eegFileIDstr = 'ptph_aeaMspmeeg_MVPA_TF';
    [feat_ph,~,~] = mvpa.extractdata(dirPreproc,condDef,eegFileIDstr,trMethod);
    save(fullfile(dirAnalysis,destFname),'feat_ph','-append');
    % feat_f and feat_ph will be 4 dimensional: 
    % 1. number of channels
    % 2. number of frequencies
    % 3. number of time samples
    % 4. number of examples
    
end


end


function [info,infoAvg,featAvg] = averagetrials(info,feat,condDef,varargin)
% Averages trials of the same condition respecting original session and day
% assginment

p = inputParser;

addRequired(p,'info');
addRequired(p,'feat');
addRequired(p,'condDef');
addParameter(p,'reuseInfo',false);

parse(p,info,feat,condDef,varargin{:});

info = p.Results.info;
feat = p.Results.feat;
condDef = p.Results.condDef;
reuseInfo = p.Results.reuseInfo;


% Assigning days
nDays = 3;
nSessionsPerDay = 20;
% Sessions
% sessions = unique(info.session);
sessions = (1:nDays*nSessionsPerDay)';
% dayIdx = repmat(1:nDays,size(sessions,1)/nDay,1);
dayIdx = repmat(1:nDays,nSessionsPerDay,1);
dayIdx = dayIdx(:);
info.day = info.session;
for i = 1:nDays, info.day(ismember(info.session,sessions(dayIdx == i))) = i; end
% Finding condition labels
AVconds = condDef.condition(~isnan(condDef.locationAuditory) & ~isnan(condDef.locationVisual));
Aconds = condDef.condition(isnan(condDef.locationVisual));
Vconds = condDef.condition(isnan(condDef.locationAuditory));
% Assigning session type variable
info.sessionType = cell(size(info.day));
info.sessionType(ismember(info.condition,Aconds)) = repmat({'A'},size(find(ismember(info.condition,Aconds))));
info.sessionType(ismember(info.condition,Vconds)) = repmat({'V'},size(find(ismember(info.condition,Vconds))));
info.sessionType(ismember(info.condition,AVconds) & info.task == 1) = repmat({'AV_a'},size(find(ismember(info.condition,AVconds) & info.task == 1)));
info.sessionType(ismember(info.condition,AVconds) & info.task == 2) = repmat({'AV_v'},size(find(ismember(info.condition,AVconds) & info.task == 2)));
info.sessionType = categorical(info.sessionType);
% Assigning new session labels
% The table must be sorted by session, condition and hand in order to keep 
% roughly the same amount of trials per condition in the averaged trials
% and consistently average both left and right response hand trials. 
nNewSessionPerDay = 4;
sessionTypes = unique(info.sessionType);
info.regroupInput = [info.session,info.iTrialSession];
    % Function for regrouping the trials for averaging
    function out = regroup(x)
        out = x;
        x = x(:,1);
        s = unique(x);
        ss = size(s,1);
        if  ss > nNewSessionPerDay
            tt = x;
            idx = zeros(ss,1);
            for ii = 1:ss, idx(ii) = mod(ii,nNewSessionPerDay); end
            idx(idx == 0) = nNewSessionPerDay;
            for ii = 1:ss, tt(ismember(x,s(ii,:))) = idx(ii); end
            out = [out,tt];
        else
            % In this case each new group is formed using trials from all
            % sessions. 
            tt = mod((1:size(x,1))',nNewSessionPerDay);
            out = [out,tt+1];
        end
        out = {out};
    end
if ~reuseInfo
    info = sortrows(info,{'session','condition','hand'}); 
    temp = varfun(@regroup,info,'InputVariables','regroupInput','GroupingVariables',...
                  {'day','sessionType'});
    sessionsRegrouped = cell2mat(temp.Fun_regroupInput);
    info.regroupInput = [];
    % Making sure the order of variables match in the table and array
    info = sortrows(info,{'day','sessionType','session','condition','hand'}); 
    temp = table2array(info(:,{'session','iTrialSession'}));
    if ~isequal(temp,sessionsRegrouped(:,1:2))
        error('mvpa:prepdata:averagetrials:exampleOrderMismatch',...
              'The order of examples doesn''t match! ');
    end
    info.newSession = sessionsRegrouped(:,3);
    info = sortrows(info,{'session','iTrialSession'});
    % Making new session values different for each day
    for i = 1:nDays
        info.newSession(info.day == i) = info.newSession(info.day == i)+(nNewSessionPerDay*(i-1));
    end
else
    if ~ismember('newSession',info.Properties.VariableNames)
        error('mvpa:prepdata:averagetrials:missingVariable',...
              'newSession must be present in info if reuseInfo is true! ');
    end
end


% Selecting conditions which are actually present
condsPresent = unique(info.condition);
condDef = condDef(ismember(condDef.condition,condsPresent),:);
% Averaging trials for each new session and condition
nConds = size(condDef,1);
featSize = num2cell(size(feat));
dimToAvg = numel(featSize);
featAvg = NaN(featSize{1:dimToAvg-1},nNewSessionPerDay*nDays*nConds);
for i = 1:(nNewSessionPerDay*nDays)
    temp = NaN(featSize{1:dimToAvg-1},nConds);
    for j = 1:nConds
        actCond = condDef.condition(j);
        isExample = info.newSession == i & info.condition == actCond;
        if all(~isExample)
            warning('mvpa:prepdata:missingExample',...
                'No example was found for session %d condition %d!',i,j);
        end
        if dimToAvg == 3
            temp(:,:,j) = mean(feat(:,:,isExample),dimToAvg);
        elseif dimToAvg == 4
            temp(:,:,:,j) = mean(feat(:,:,:,isExample),dimToAvg);
        end
    end
    startIdx = 1+(nConds*(i-1));
    endIdx = nConds+(nConds*(i-1));
    if dimToAvg == 3
        featAvg(:,:,startIdx:endIdx) = temp;
    elseif dimToAvg == 4
        featAvg(:,:,:,startIdx:endIdx) = temp;
    end
end
% Creating info for the averaged trials
infoAvg = repmat(condDef,nNewSessionPerDay*nDays,1);
infoAvg.Properties.VariableNames = {'condition','task','locA','locV','relV'};
temp = repmat(1:nNewSessionPerDay*nDays,nConds,1);
infoAvg.session = temp(:);
infoAvg.resp = NaN(size(infoAvg.session));
infoAvg = infoAvg(:,{'session','condition','task','locA','locV','relV','resp'});
% Checking if the respective sizes match
if size(infoAvg,1) ~= size(featAvg,dimToAvg)
    error('mvpa:prepdata:averagetrials:exampleNumberMismatch',...
        'Inconsistent number of examples in info and feat! ');
end

end


function tOut = assignnewsessions(info,nNewSession)

tIn = unique(info(:,{'day','session','task'}));
nDay = max(tIn.day);
nNewSessionPerDay = nNewSession/nDay;
newSession = NaN(size(tIn,1),1);
nSessionInDay = sum(tIn.day == 1);

for iDay = 1:nDay
    
    tempNewSession = [];
    for i = 1:nSessionInDay/nNewSessionPerDay
        tempNewSession = [tempNewSession;randperm(nNewSessionPerDay)']; %#ok<AGROW>
    end
    
    tempNewSession = tempNewSession+((iDay-1)*nNewSessionPerDay);
    startIdx = ((iDay-1)*nSessionInDay)+1;
    endIdx = ((iDay-1)*nSessionInDay)+nSessionInDay;
    
    newSession(startIdx:endIdx) = tempNewSession;
    
end

% Making sure, that for a given new session equal number of sessions are
% going to be assigned from both task types. 
tOut = sortrows(tIn,{'day','task'});
tOut.newSession = newSession;
tOut = sortrows(tOut,{'session'});

end
