function varargout = run_mvpa_fun_group(fun,trMethod,matchStr,varargin)
% Runs mvpa functions on the group of subjects

%% Parsing input, checking matlab
p = inputParser;

validTrMethods = {'sample-wise-sm','sample-wise-sm-avg','sample-wise-bp',...
    'sample-wise-sm-corr1','sample-wise-sm-corr2','sample-wise-sm-avg-corr2',...
    'sample-wise-bp-avg','sample-wise-source','sample-wise-source-avg'};
validExpStages = {'pilot_1','pilot_2','pilot_3','pilot_4','final'};

addRequired(p,'fun',@(x) validateattributes(x,{'function_handle'},{'scalar'}));
addRequired(p,'trMethod',@(x) any(validatestring(x, ...
                                                 validTrMethods)));
addRequired(p,'matchStr',@ischar);
addParameter(p,'expStage','final',@(x) any(validatestring(x,validExpStages)));
addParameter(p,'funInput',{},@(x) validateattributes(x,{'cell'},{'vector','row'}));
addParameter(p,'nFilesExpected',[],@(x) validateattributes(x,{'numeric'},{'scalar'}));
addParameter(p,'subFolder',[],@(x) validateattributes(x,{'char'},{'nonempty'}));

parse(p,fun,trMethod,matchStr,varargin{:});

fun = p.Results.fun;
trMethod = p.Results.trMethod;
expStage = p.Results.expStage;
fileMatchStr = p.Results.matchStr;
funInput = p.Results.funInput;
nFilesExpected = p.Results.nFilesExpected;
subFolder = p.Results.subFolder;

% Finding subject folders
saveDf = cd(DEC_2_setupdir(expStage,'anal_eeg'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

% Collecting subject level mvpares datasets
dirID = 'anal_eeg_sub_mvpa';
miscInput = {};
% Preparations specific to certain functions
if strcmp(func2str(fun),'mvpa.addAVmodelCorrelations')
    % Expected Number of files per subject
    if isempty(nFilesExpected)
        nFilesExpected = 1;
    end
    filePathList = collectFiles(subjList,expStage,dirID,fileMatchStr,trMethod,nFilesExpected,subFolder);
    dirID = 'anal_behav_sub';
    fileMatchStrBehav = 'BEHAV_ANAL_.*';
    trMethod = '';
    pathFmriFile = fullfile(DEC_2_setupdir(expStage,'study_root'),'docs',...
        'Tim_model_free','wav_tim_group.mat');
    behavFilePathList = collectFiles(subjList,expStage,dirID,fileMatchStrBehav,trMethod,nFilesExpected,subFolder);
    if isempty(behavFilePathList)
        miscInput = [repmat({'pathFmriFile'},size(filePathList)),...
            repmat({pathFmriFile},size(filePathList))];
    else
        miscInput = [repmat({'pathFmriFile'},size(filePathList)),...
            repmat({pathFmriFile},size(filePathList)),...
            repmat({'pathBehavFile'},size(filePathList)),behavFilePathList];
    end
elseif strcmp(func2str(fun),'mvpa.mergeMvpaRes')
    % Expected Number of files per subject
    if isempty(nFilesExpected)
        nFilesExpected = 2;
    end
    filePathList = collectFiles(subjList,expStage,dirID,fileMatchStr,trMethod,nFilesExpected,subFolder);
else
    % Expected Number of files per subject
    if isempty(nFilesExpected)
        nFilesExpected = 2;
    end
    filePathList = collectFiles(subjList,expStage,dirID,fileMatchStr,trMethod,nFilesExpected,subFolder);
end

% Executing function on all subject level data
if strcmp(func2str(fun),'mvpa.averageMvpaRes')
    matchStr = '([cenr]{2}_)(tr[A-Za-z-() ]+_)([a-zA-Z]+_)(?:[0-9-]+_)(roi-[a-zA-Z0-9-]+_)?(gen[A-Za-z-() ]+)';
    id = regexp(filePathList{1},matchStr,'tokens');
    id = id{:};
    id = cat(2,id{:});
    % Getting rid of special caracters to create valid file name
    id = regexprep(id, '[/\*:?"<>|]*',' ');
    avgFileName = ['gr_',id,'.mat'];
    
    I.pathAveragedFile = fullfile(DEC_2_setupdir(expStage,'anal_eeg_group_mvpa'),...
                                  trMethod,avgFileName);
    if ~exist(fullfile(DEC_2_setupdir(expStage, ...
                                      'anal_eeg_group_mvpa'),trMethod),'file')
        mkdir(fullfile(DEC_2_setupdir(expStage,'anal_eeg_group_mvpa'),trMethod));
    end
    I.pathFilesToAverage = filePathList;
    varargout{1} = mvpa.averageMvpaRes(I);
elseif strcmp(func2str(fun),'mvpa.mergeMvpaRes')
    matchStr = '([cenr]{2}_)(tr[A-Za-z-() ]+_)([a-zA-Z]+_)(?:[0-9-]+_)(roi-[a-zA-Z0-9-]+_)?(gen[A-Za-z-() ]+)';
    id = regexp(fileMatchStr,matchStr,'tokens');
    id = id{:};
    id = cat(2,id{:});
    % Getting rid of special caracters to create valid file name
    id = regexprep(id, '[/\*:?"<>|]*',' ');
    for i = 1:numel(subjList)
        
        subID = subjList{i};
        
        mergedFileName = [id,'.mat'];
        
        I.pathMergedFile = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID),...
                                    trMethod,mergedFileName);
        if ~exist(fullfile(DEC_2_setupdir(expStage, ...
                                          'anal_eeg_group_mvpa'),trMethod),'file')
            mkdir(fullfile(DEC_2_setupdir(expStage,'anal_eeg_group_mvpa'),trMethod));
        end
        I.pathFilesToMerge = filePathList(i,:);
        varargout{i} = mvpa.mergeMvpaRes(I); %#ok
        
    end
    
else
    for i = 1:size(filePathList,1)
        temp = mvpares(filePathList{i});
        if ~isempty(miscInput)
            % This cell array must be a row vector
            args = cat(2,funInput,miscInput(i,:));
        else
            args = funInput;
        end
        varargout{i} = fun(temp,args{:}); %#ok
    end
end

end

function filePathList = collectFiles(subjList,expStage,dirID,fileMatchStr,trMethod,nFilesExpected,subFolder)

filePathList = {};
index = 1;
for i = 1:size(subjList,1)
    saveDf = cd(fullfile(DEC_2_setupdir(expStage,dirID,subjList{i}),trMethod,subFolder));
    fileList = cellstr(ls);
    matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
    if sum(matchID) == 0
        warning('No file, skipping subject %s! ',subjList{i});
        cd(saveDf);
        continue;
    elseif sum(matchID) > nFilesExpected
        warning('More files than needed, skipping subject %s! ',subjList{i});
        cd(saveDf);
        continue;
    else
        fileName = fileList(matchID);
        if iscolumn(fileName)
            fileName = fileName';
        end
    end
    temp = cellfun(@fullfile,repmat({pwd},size(fileName)),fileName,'UniformOutput',false);
    filePathList = cat(1,filePathList,temp);
    index = index + 1;
    cd(saveDf);
end

if isrow(filePathList)
    filePathList = filePathList';
end

end
