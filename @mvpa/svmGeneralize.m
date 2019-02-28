function [mvparesObj,filePath] = svmGeneralize(I)
% Method for generalizing a support vector machine model on timeseries data
% 
% USAGE:
%   [mvparesObj,filePath] = svmGeneralize(I)
% INPUT: 
%   I (struct): input settings. 
%       Required fields:
%           dir_dataFile (string): path to the MVPA dataset
%           dir_trainFile (string): path to the trained mvpares dataset
%           gen_cond (string): generalization condition
%           gen_time (string): indicating the generalization time. Possible 
%               values: 'tr','tr_x_tr' for trainin time and training time 
%               by training time respectively
%       Optional fields:
%           customTag (string): arbitrary tag to be added to the name of
%               the output file
%           progrMonitor (logical): wether to display a progress monitor,
%              default: true
% OUTPUT:
%   filePath: full path to the output mvpa result structure saved on disc 
%       as a .mat file. The file contains the variables/fields of the
%       training mvpa result file and the following fields with the 
%       generalization results:
%       gen_examples: table of information about the examples used for 
%           generalization. Number of rows = number of examples, 
%       gen_accuracies: M x N x O x P array with 
%           M = 3 (accuracy,MSE,R2), 
%           N = number of models per sample
%           O = number of training samples
%           P = number of generalization samples per training sample
%       gen_predlabels: M x N x O x P array with 
%           M = maximum number of generalization examples during the 
%               cross-validation padded with NaNs if the actual number of 
%               test examples is smaller than that)
%           N = number of models per sample
%           O = number of samples 
%           P = number of generalization samples per training sample
%       The info variable is extended with further generalization details.
%   mvparesObj: The above results loaded as an mvpares object. 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

%% Parsing input, checking matlab
p = inputParser;

% Valid values for certain inputs.
validGenConds = {'A','V','AV','AxorV','pseudoAV','c','ci','i','a','v','av','hrel','lrel',...
    'lh','rh','la','ha'};
validGenTimes = {'tr','tr_x_tr'};
validTrGenCondMatch = {'no','partial','full'};

% Defining inputs.
addParameter(p,'customTag','',@ischar);
addParameter(p,'dir_trainFile','',@(x)exist(x,'file'));
addParameter(p,'dir_dataFile','',@(x)exist(x,'file'));
addParameter(p,'gen_cond','',@(x)all(ismember(regexp(x,'-','split'),validGenConds)));
addParameter(p,'gen_time','',@(x)any(validatestring(x,validGenTimes)));
addParameter(p,'tr_genCondMatch','',@(x)any(validatestring(x,validTrGenCondMatch)));
addParameter(p,'progrMonitor',true,@islogical);

% Parsing inputs.
parse(p,I);

% Assigning inputs to variables.
customTag = p.Results.customTag;
dirTrainFile = p.Results.dir_trainFile;
dirDataFile = p.Results.dir_dataFile;
genCond = p.Results.gen_cond;
genTimeStr = p.Results.gen_time;
tr_genCondMatch = p.Results.tr_genCondMatch;
progrMonitor = p.Results.progrMonitor;

% Custom checkings after parsing
requiredVars = {'dir_trainFile','dir_dataFile','gen_cond',...
                'gen_time'};
if any(ismember(requiredVars,p.UsingDefaults))
    error('mvpa:svmGeneralize:missingInput',...
        'All required parameters must be specified!');
end

%% Checking... 
% matlab version
if verLessThan('matlab', '8.3.0.532')
    error('mvpa:svmGeneralize:notSupportedMatlabVersion',...
        'MATLAB 2014a or newer is required! ');
end

% parallel pool 
currPool = gcp('nocreate');
if isempty(currPool)
    error('mvpa:svmGeneralize:missingParallelPool',...
        'Please open a parallel pool in order to run the function in parallel mode!');
end

% libsvm 
if isempty(regexp(path,'libsvm','once'))
    error('mvpa:svmGeneralize:missingToolbox','libsvm not found in path!');
end 

% parfor_progress
if progrMonitor
    if isempty(regexp(path,'parfor_progress','once'))
        error('mvpa:svmGeneralize:missingToolbox','parfor_progress not found in path!');
    end
end

%% Starting the timer
cStartFun = clock;
ds = datestr(now);
fprintf('\nsvmGeneralize%s%s\n%s\n',repmat(' ',1,67-length(ds)),ds,repmat('-',1,80));

%% Loading and checking data files...
trainM = matfile(dirTrainFile);
dataM = matfile(dirDataFile);

%% Defining importand variables
trInfo = trainM.info;
cvScheme = trInfo.cv_scheme;
scMethod = trInfo.sc_method;
scParam = trainM.tr_scParam;
[~,name,ext] = fileparts(dirDataFile);
genDataFile = [name,ext];
if isfield(trInfo,'tr_data_file')
    trDataFile = trInfo.tr_data_file;
else
    % If the training data file is unknown we assume it is the same as the
    % generalization as this is the most likely case and also this is for
    % backward compatibility
    trDataFile = genDataFile;
end
trIsExample = trInfo.tr_isExample;
trMethod = trInfo.tr_method;
trTimePoints = trInfo.tr_timePoints;
trCond = trInfo.tr_cond;
trLabel = trInfo.tr_label;
trModels = trainM.tr_models;
if ~isempty(regexp(trCond,'pseudoAV','once'))
    trGroupings = trainM.tr_groupings_AxorV;
else
    trGroupings = trainM.tr_groupings;
end
% Number of training time points
nTrTimePoints = size(trTimePoints,2);
% Determining generalization label
if any(ismember(trLabel,{'finger','resp','hand','disp',...
                        'relV','task','LvsRa','LvsRv'}))
    genLabelInput = trLabel;
elseif strcmp(trLabel,'sA_hat_indep')
    genLabelInput = 'sAhatInd';
elseif strcmp(trLabel,'sV_hat_indep')
    genLabelInput = 'sVhatInd';
elseif strcmp(trLabel,'s_hat_comm')
    genLabelInput = 'sHatComm';
elseif strcmp(trLabel,'sA_hat_and_sV_hat')
    genLabelInput = 'sAhatANDsVhat';
else
    genLabelInput = 'stim';
end

% Cross-validation specific variables
if strcmp(cvScheme,'kf')
    nKFrep = trInfo.nKFrep;
    k = trInfo.k;
else
    nKFrep = [];
    k = [];
end

% Dataset to generalize to
condDef = dataM.condDef;
dataMisc = dataM.misc;
if any(strcmp(trMethod,{'sample-wise-beta','sample-wise-beta-rl','sample-wise-sm-beta'}))
    if strcmp(genLabelInput,'stim')
        infoStr = 'info_stim';
    elseif strcmp(genLabelInput,'resp')
        infoStr = 'info_resp';
    elseif strcmp(genLabelInput,'hand')
        error('mvpa:svmGeneralize:missingFunctionality',...
            'This analysis is not yet implemented!');
    end
else
    infoStr = 'info';
end
dataInfo = dataM.(infoStr);

%% Choosing the subset of data which is used for generalizing.
if any(strcmp(trMethod,{'sample-wise-beta','sample-wise-beta-rl','sample-wise-sm-beta'})) && strcmp(genLabelInput,'resp')
    
    % In this case we don't really have stimulus conditions as all of them
    % are represented in the glm, so we set the trCond manually and use all
    % the examples available.
    if ~strcmp(genCond,'AV-ci-av')
        warning('mvpa:svmGeneralize:automaticOptionSetting',...
            'Changing training condition automatically to AV-ci-av!');
        genCond = 'AV-ci-av';
    end
    genIsExample = true(size(dataInfo,1),1);
    genLabel = 'resp';
    
else
    
    genIsExample = mvpa.selectexamples(genCond,condDef,dataInfo);
    if sum(genIsExample) == 0
        error('mvpa:svmGeneralize:requestedDataNotPresent',...
            'No examples of the specified condition were found in the dataset!');
    end
    genLabel = mvpa.selectlabel(genLabelInput,genCond);
    
end

% Determining if the training condition and the generalization condition
% are identical. 
if ~strcmp(trDataFile,genDataFile)
    % If the training and generalization data files are different
    % it is possible to specify exactly the relationship between
    % the training and generalization examples. By default this is
    % 'no'. If other value than 'no' is specified, the training and
    % generalization datasets must have identical rows (only true
    % for sample-wise-sm-avg and sample-wise-sm-avg-corr data). 
    if isempty(tr_genCondMatch)
        tr_genCondMatch = 'no';
    end
elseif all(~ismember(find(genIsExample),find(trIsExample)))
    % If the training and generalization example sets are disjoint
    tr_genCondMatch = 'no';
elseif any(~ismember(find(trIsExample),find(genIsExample)))
    % If the training and generalization set is not disjoint but not all
    % training examples are present in the generalization example set. 
    error('mvpa:svmGeneralize:datasetMismatch',...
        ['All training examples should be present among the generalization ',...
        'examples if they are not completely different!']);
elseif all(ismember(find(genIsExample),find(trIsExample)))
    % If the generalization example set matches the training example set
    % completely
    tr_genCondMatch = 'full';
else
    % If the generalization example set contains examples which are not
    % part of the training example set
    tr_genCondMatch = 'partial';
end


%% Extracting the appropriate set of examples and features for generalizing
% genExamples_feat is a 3D matrix with 
% size(1) = nFeatures
% size(2) = nTimePoints
% size(3) = nExamples
% Feature scaling is also done here. 
genTimePoints = trTimePoints;
if strcmp(trMethod,'sample-wise-tf')
    genExamples_feat = mvpa.extractfeatures(dataM,genIsExample,trMethod,genLabelInput,genTimePoints,...
        'tfFreq',trInfo.tf_freq,'tfType',trInfo.tf_type,'progrMonitor',progrMonitor);
elseif ismember(trMethod,{'sample-wise-bp','sample-wise-bp-avg'})
    genExamples_feat = mvpa.extractfeatures(dataM,genIsExample,trMethod,genLabelInput,genTimePoints,...
        'bpFreq',trInfo.bp_freq,'progrMonitor',progrMonitor);
elseif ismember(trMethod,{'sample-wise-source','sample-wise-source-avg'})
    genExamples_feat = mvpa.extractfeatures(dataM,genIsExample,trMethod,genLabelInput,genTimePoints,...
        'trROI',trInfo.tr_roi,'progrMonitor',progrMonitor);
else
    genExamples_feat = mvpa.extractfeatures(dataM,genIsExample,trMethod,genLabelInput,genTimePoints,...
        'progrMonitor',progrMonitor);
end
genExamples_info = dataInfo(genIsExample,:);

% Determining the number of generalization time points. 
if strcmp(genTimeStr,'tr')
    nGenTimePoints = 1;
elseif strcmp(genTimeStr,'tr_x_tr')
    nGenTimePoints = nTrTimePoints;
end

% The maximum number of predicted labels over all models. 
switch tr_genCondMatch
    case 'full'
        % The grouping is the same as for training
        genGroupings = trGroupings;
    case 'partial'
        % Finding examples not present in the training example set
        genIsExampleTrNotExample = genIsExample & ~trIsExample;
        % Assigning groups to the above examples according to the cv scheme
        if strcmp(cvScheme,'kf')
            misc.k = k;
            misc.nKFrep = nKFrep;
        end
        misc.nTimePoints = nTrTimePoints;
        tempGroupings = mvpa.assigngroupings(cvScheme,...
            dataInfo(genIsExampleTrNotExample,:),misc);
        % Merging the grouping with the training grouping
        genGroupings = NaN(size(dataInfo,1),size(tempGroupings,2),...
            size(tempGroupings,3));
        genGroupings(trIsExample,:,:) = trGroupings;
        genGroupings(genIsExampleTrNotExample,:,:) = tempGroupings;
        genGroupings = genGroupings(genIsExample,:,:);
    case 'no'
        % Assigning groups to the new examples according to the cv scheme
        % Assigning groups to the above examples according to the cv scheme
        if strcmp(cvScheme,'kf')
            misc.k = k;
            misc.nKFrep = nKFrep;
        end
        misc.nTimePoints = nTrTimePoints;
        genGroupings = mvpa.assigngroupings(cvScheme,genExamples_info,misc);
end

% We need to generate the pseudoAV trials from A and V trials here
% and this will affect the number of examples in the dataset
if ~isempty(regexp(genCond,'pseudoAV','once'))
    [genExamples_info,genExamples_feat,genGroupings] = ...
        mvpa.generatePseudoAvData(genExamples_info,genExamples_feat,genGroupings,genCond,condDef);
end

% Selecting generalization labels
if ~ismember(genLabel,genExamples_info.Properties.VariableNames)
    error('mvpa:svmGeneralize:invalidVariableValue',...
        'The specified generalization label is not part of the dataset!');
else
    genExamples_lab = genExamples_info.(genLabel);
end
if any(isnan(genExamples_lab))
    error('mvpa:svmGeneralize:invalidVariableValue',...
        'Generalization labels can''t be NaNs! ');
end
sizeGenExamples_feat = size(genExamples_feat);

% The maximum number of test examples in the cross-validation
% scheme.
nPredLabels = maxnival(genGroupings(:,1,1));

% Number of models per sample and in total. 
nModelsPerTimePoint = size(trModels,3);

% Pre-allocating arrays for data collection. 
gen_predlabels = NaN(nPredLabels,nModelsPerTimePoint,nTrTimePoints,nGenTimePoints);
gen_accuracies = NaN(3,nModelsPerTimePoint,nTrTimePoints,nGenTimePoints);

% Initializing progress monitor
cStartGen = clock;
fprintf('Generalizing...\n');
if progrMonitor
    parfor_progress(nTrTimePoints*nGenTimePoints);
end

% Iterating through training samples
parfor iTrTimePoint = 1:nTrTimePoints
    
    % Selecting the appropriate sample of models and scalin parameters. 
    actTrTimePointModels = trModels(:,:,:,iTrTimePoint);
    actTrTimePointScParams = scParam(:,:,:,iTrTimePoint);
    
    % Pre allocating arrays for collecting data for a given sample. 
    actTrTimePointGenPredlab = NaN(nPredLabels,nModelsPerTimePoint,1,nGenTimePoints);
    actTrTimePointGenAcc = NaN(3,nModelsPerTimePoint,1,nGenTimePoints);
    
    for iGenTimePoint = 1:nGenTimePoints
        
        % Extracting the actual sample. This yields a matrix with
        % number of rows = number of examples,
        % number of columns = number of features.
        if strcmp(genTimeStr,'tr')
            actGenTimePointFeats = squeeze(genExamples_feat(:,iTrTimePoint,:))';
        else
            actGenTimePointFeats = squeeze(genExamples_feat(:,iGenTimePoint,:))'; 
        end
        
        if any(size(actGenTimePointFeats) ~= [sizeGenExamples_feat(3),sizeGenExamples_feat(1)])
            error('mvpa:svmGeneralize:datasetMismatch','Incorrect feature matrix sizes.');
        end
        
        % We do the generalization according to the cross-validation
        % scheme.
        
        % Preallocating arrays for test labels and features
        teLabs = NaN(nPredLabels,nModelsPerTimePoint);
        teFeats = NaN(nPredLabels,size(actGenTimePointFeats,2),nModelsPerTimePoint);
        
        if strcmp(cvScheme,'kf')
            for iKFrep = 1:nKFrep
                actGrouping = genGroupings(:,iKFrep,iTrTimePoint);
                for iFold = 1:k
                    iModel = ((iKFrep-1)*k)+iFold;
                    actFoldExamples = actGrouping == iFold;
                    actNumExamples = sum(actFoldExamples);
                    teLabs(1:actNumExamples,iModel) = genExamples_lab(actFoldExamples);
                    teFeats(1:actNumExamples,:,iModel) = ...
                        mvpa.scalefeatures(actGenTimePointFeats(actFoldExamples,:),...
                        scMethod,actTrTimePointScParams(:,:,iModel));
                end
            end
        elseif strcmp(cvScheme,'loso')
            actGrouping = genGroupings(:,iTrTimePoint);
            for iModel = 1:nModelsPerTimePoint
                actModelExamples = actGrouping == iModel;
                actNumExamples = sum(actModelExamples);
                teLabs(1:actNumExamples,iModel) = genExamples_lab(actModelExamples);
                teFeats(1:actNumExamples,:,iModel) = ...
                    mvpa.scalefeatures(actGenTimePointFeats(actModelExamples,:),...
                    scMethod,actTrTimePointScParams(:,:,iModel));
            end
        end
        
        for iModel = 1:nModelsPerTimePoint
            
            actTeLabs = teLabs(:,iModel);
            actTeLabs = actTeLabs(~isnan(actTeLabs));
            actTeFeats = teFeats(:,:,iModel);
            actTeFeats = actTeFeats(~isnan(actTeLabs),:);
            
            [plabs,acc,~] = svmpredict(actTeLabs,actTeFeats,...
                mvpa.mat2mdl(actTrTimePointModels(:,:,iModel)));
            
            % Padding with NaNs if necesary
            if size(plabs,1) < nPredLabels
                plabs = [plabs;NaN(nPredLabels-size(plabs,1),1)];
            end
            
            actTrTimePointGenPredlab(:,iModel,1,iGenTimePoint) = plabs;
            actTrTimePointGenAcc(:,iModel,1,iGenTimePoint) = acc;
            
        end
        
        % Advancing Progress monitor
        if progrMonitor
            parfor_progress;
        end
    end
    
    gen_predlabels(:,:,iTrTimePoint,:) = actTrTimePointGenPredlab;
    gen_accuracies(:,:,iTrTimePoint,:) = actTrTimePointGenAcc;
    
end

% Finalizing progress monitor.
if progrMonitor
    parfor_progress(0);
end
fprintf('Generalization elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(clock,cStartGen)/86400,'dd HH:MM:SS'));

% If the generalization condition does not have an unambiguous label we do
% not save generalization label and accuracies. 
genCondsCell = regexp(genCond,'-','split');
if strcmp(genLabelInput,'stim') && (all(ismember({'ci','av'},genCondsCell)) || ...
        all(ismember({'i','av'},genCondsCell)) || all(ismember({'AV','ci'},genCondsCell)))
    gen_accuracies = NaN;
    genLabel = '';
end

% Saving the output into the struct. 
% Updating the info field
trInfo.gen_cond = genCond;
trInfo.gen_elapsedTime = etime(clock,cStartFun);
trInfo.gen_isExample = genIsExample;
trInfo.gen_label = genLabel;
trInfo.gen_time = genTimeStr;
[~,name,ext] = fileparts(dirDataFile);
trInfo.gen_data_file = [name,ext];
trInfo.gen_data_fs = dataMisc.fs;
% Saving computed data
genM = struct();
genM.gen_accuracies = gen_accuracies;
genM.gen_examples = genExamples_info;
genM.gen_groupings = genGroupings;
genM.gen_predlabels = gen_predlabels;
genM.tr_accuracies = trainM.tr_accuracies;
genM.tr_examples = trainM.tr_examples;

% Closing the training file. 
trainM = [];

% Ordering the sturcture fields. 
trInfo = orderfields(trInfo);
genM.info = trInfo; %#ok<*STRNU>

% Substring for generalization label
genLabelStr = genLabelInput;

% Generating output file name. 
if strcmp(genTimeStr,'tr')
    genTimeStr = 't';
elseif strcmp(genTimeStr,'tr_x_tr')
    genTimeStr = 'txt';
end

dateString = datestr(now,'yymmddHHMMSS');
if ~isempty(customTag), customTag = [customTag,'-']; end
fileNameAppend = sprintf('gen-%s%s_%s_%s_%s',customTag,genCond,genLabelStr,genTimeStr,...
    dateString);
extStart = strfind(dirTrainFile,'.mat');
filePath = [dirTrainFile(1:extStart-(length(dateString)+1)),...
    fileNameAppend,dirTrainFile(extStart:end)];

% Saving the outputs. 
save(filePath,'-struct','genM','-v7.3');
mvparesObj = mvpares(filePath);

end