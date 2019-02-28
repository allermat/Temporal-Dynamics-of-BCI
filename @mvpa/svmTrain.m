function [mvparesObj,filePath] = svmTrain(I)
% Method for training a support vector machine model on timeseries data
%
% USAGE:
%   [mvparesObj,filePath] = svmTrain(I)
% INPUT: 
%   I (structure): input settings with the following fields:  
%       cv_scheme: crossvalidation scheme, valid schemes: 'kf','loso','loxo'
%       dir_analysis: full path of the analysis directory. 
%       dir_dataFile: full path of the data file. 
%       dir_utils: full path of the utilities directory. 
%       k: number of folds in the k-fold crossvalidation 
%       nKFrep: number of repetitions of the k-fold crossvalidation
%       params: 1x2 cell array containing the vectors of possible parameters to
%           check during the gridsearch. 
%       subID:
%       svm_type: cc - C-SVC, nc - nu-SVC, er - epsilon-SVR, nr - nu-SVR
%       tr_cond: hyphen separated list of training conditions, valid
%           conditions: 'A','V','AV','AxorV','c','ci','i','a','v','av','lh','rh'
%       tr_label: the label for training, valid labels: 'hand','resp',
%           'cond','sAhatInd','sVhatInd','sHatComm','sAhatANDsVhat',
%           'LvsRa','LvsRv'
%       tr_method: training method, valid methods: 
%           'sample-wise','sample-wise-beta','sample-wise-beta-rl',...
%           'sample-wise-bp','sample-wise-bp-avg','sample-wise-rl',...
%           'sample-wise-sc','sample-wise-mult','sample-wise-sm',...
%           'sample-wise-tf'
%       tr_timePoints: cell array of samples in ms to be used. 
%           If there are multiple time points specified in one cell, than the
%           features belonging to these time points are going to be treated as 
%           one set. 
%       progrMonitor (logical): wether to display a progress monitor 
% OUTPUT:  
%   filePath: full path to the output mvpa result structure saved on disc 
%       as a .mat file. The file contains the variables/fields below. 
%       The dimensions of the variables varies depending on the 
%       cross-validation method. 
%       
%       K-fold CV:
%       ==========
%       info: 1x1 struct with various information
%           
%       tr_examples: table of information about the examples used for 
%           training. Number of rows = number of examples,  
%       tr_grids: 1x2 struct array with the following fields
%           pass, grids, param1, param2 (for SVRs), bestPar
%       tr_groupings: M x N x O array with 
%           M = number of examples, 
%           N = I.nKFrep
%           O = number of time points 
%       tr_models: M x N x O x P array with M,N determined by the maximum 
%           size of the matrices generated from the output models of 
%           svmtrain, O = number of models per time point, P = number of
%           time points.
%       tr_scParam: M x N x O x P array with M = number of features, 
%           N = 2 (mean and std for Z-scoring), O = number of models per 
%           sample, P = number of samples. 
%       
%       Leave One Session Out (LOSO) CV:
%       ================================
%       info: 1x1 struct with various information about training
%           
%       tr_examples: table of information about the examples used for 
%           training. Number of rows = number of examples,  
%       tr_grids: 1x2 struct array with the following fields: 
%           pass, grids, param1, param2 (for SVRs), bestPar
%       tr_groupings: M x N array with 
%           M = number of examples 
%           N = number of time points 
%       tr_models: M x N x O x P array with M,N determined by the maximum 
%           size of the matrices generated from the output models of 
%           svmtrain, O = number of models per time point, P = number of
%           time points.
%       tr_scParam: M x N x O x P array with M = number of features, 
%           N = 2 (mean and std for Z-scoring), O = number of models per 
%           time point, P = number of time points.
%
%   mvparesObj: The above results contained in an mvpares object. 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

%% Parsing input
p = inputParser;

% Valid values for certain inputs. 
validBpFreq = {'d','t','a','b','gl','gh'};
validCVschemes = {'kf','loso','loxo'};
validSCmethods = {'one-f','one-e','Z-f','Z-e','Z-ef','Z-fe'};
validSVMtypes = {'cc','er','nc','nr'};
validROIs = {'A','IPS','V1-3'};
validTfTypes = {'f','ph','fph'};
validTrConds = {'A','V','AV','AxorV','pseudoAV','c','ci','i','a','v','av','hrel','lrel',...
    'lh','rh','la','ha'};
validTrLabels = {'finger','hand','resp','stim','sAhatInd','sVhatInd','sHatComm',...
    'sAhatANDsVhat','LvsRa','LvsRv'};
validTrMethods = {'sample-wise','sample-wise-avg','sample-wise-beta',...
    'sample-wise-beta-rl','sample-wise-bp','sample-wise-bp-avg',...
    'sample-wise-rl','sample-wise-sc','sample-wise-mult','sample-wise-sm',...
    'sample-wise-sm-avg','sample-wise-sm-beta','sample-wise-sm-corr1', ...
    'sample-wise-sm-corr2','sample-wise-sm-avg-corr2','sample-wise-tf',...
    'sample-wise-source','sample-wise-source-avg'};

% Custom validation functions for certain inputs. 
    function checkTrTimePoints(x)
        if ~iscell(x) || ~isrow(x)
            error('mvpa:svmTrain:invalidInput','tr_timePoints must be a cell array of one row!');
        elseif any(~cellfun(@isrow,x)) || numel(unique(cellfun(@numel,x))) > 1
            error('mvpa:svmTrain:invalidInput',...
                'All elements of tr_timePoints must be row vectors of equal length!');
        elseif any(~cellfun(@ismember,num2cell(cellfun(@mean,x)),x))
            error('mvpa:svmTrain:invalidInput',...
                'All cells of tr_timePoints must contain their own mean!');
        end
    end

    function checkTfFreq(x)
        if ~iscell(x) || ~isrow(x) || numel(x) > 2
            error('mvpa:svmTrain:invalidInput',...
                'tf_freq must be a cell array of one row, maximum two elements!');
        elseif any(~cellfun(@isrow,x))
            error('mvpa:svmTrain:invalidInput',...
                'All elements of tf_freq must be row vectors!');
        end
    end

% Defining inputs.
addParameter(p,'cv_scheme','',@(x)any(validatestring(x,validCVschemes)));
addParameter(p,'dir_analysis','',@(x)exist(x,'dir'));
addParameter(p,'dir_dataFile','',@(x)exist(x,'file'));
addParameter(p,'k',[],@(x)validateattributes(x,{'numeric'},...
    {'scalar','integer','positive'}));
addParameter(p,'nKFrep',[],@(x)validateattributes(x,{'numeric'},...
    {'scalar','integer','positive'}));
addParameter(p,'params',{},@(x)validateattributes(x,{'cell'},{'nrows',1}));
addParameter(p,'progrMonitor',true,@islogical);
addParameter(p,'sc_method','',@(x)any(validatestring(x,validSCmethods)));
addParameter(p,'subID','',@ischar);
addParameter(p,'svm_type','',@(x)any(validatestring(x,validSVMtypes)));
addParameter(p,'bp_freq','',@(x)any(validatestring(x,validBpFreq)));
addParameter(p,'tf_freq',{},@checkTfFreq);
addParameter(p,'tf_type','',@(x)any(validatestring(x,validTfTypes)));
addParameter(p,'tr_roi','',@(x)any(validatestring(x,validROIs)));
addParameter(p,'tr_cond','',@(x)all(ismember(regexp(x,'-','split'),validTrConds)));
addParameter(p,'tr_label','',@(x)any(validatestring(x,validTrLabels)));
addParameter(p,'tr_method','',@(x)any(validatestring(x,validTrMethods)));
addParameter(p,'tr_timePoints',[],@checkTrTimePoints);

% Parsing inputs.
parse(p,I);

% Assigning inputs to variables. 
cvScheme = p.Results.cv_scheme;
dirAnalysis = p.Results.dir_analysis;
dirDataFile = p.Results.dir_dataFile;
k = p.Results.k;
nKFrep = p.Results.nKFrep;
params = p.Results.params;
progrMonitor = p.Results.progrMonitor;
scMethod = p.Results.sc_method;
subID = p.Results.subID;
svmType = p.Results.svm_type;
bpFreq = p.Results.bp_freq;
tfFreq = p.Results.tf_freq;
tfType = p.Results.tf_type;
trROI = p.Results.tr_roi;
trCond = p.Results.tr_cond;
trLabelInput = p.Results.tr_label;
trMethod = p.Results.tr_method;
trTimePoints = p.Results.tr_timePoints;

% Custom checkings after parsing
requiredVars = {'cv_scheme','dir_analysis','dir_dataFile','sc_method',...
    'subID','svm_type','tr_cond','tr_label','tr_method','tr_timePoints'};
if any(ismember(requiredVars,p.UsingDefaults))
    error('mvpa:svmTrain:missingInput',...
        'All required parameters must be specified!');
end

requiredVarsKF = {'k','nKFrep'};
if strcmp(cvScheme,'kf') && any(ismember(requiredVarsKF,p.UsingDefaults))
    error('mvpa:svmTrain:missingInput','''k'' and ''nKFrep'' must be specified!');
end

requiredVarsTF = {'tf_freq','tf_type'};
if strcmp(trMethod,'sample-wise-tf')
    if any(ismember(requiredVarsTF,p.UsingDefaults))
        error('mvpa:svmTrain:missingInput','tf_freq and tr_type must be specified!');
    elseif strcmp(tfType,'fph') && size(tfFreq,2) < 2
        error('mvpa:svmTrain:invalidInput',...
            'tf_freq must be a 1 x 2 cell array if tf_type is ''fph''!');
    end
end

requiredVarsBP = {'bp_freq'};
if strcmp(trMethod,'sample-wise-bp') && any(ismember(requiredVarsBP,p.UsingDefaults))
    error('mvpa:svmTrain:missingInput','''bp_freq'' must be specified!');
end

requiredVarsBP = {'tr_roi'};
if ismember(trMethod,{'sample-wise-source','sample-wise-source-avg'}) ...
        && any(ismember(requiredVarsBP,p.UsingDefaults))
    error('mvpa:svmTrain:missingInput','''tr_roi'' must be specified!');
end

trCondsCell = regexp(trCond,'-','split');
if strcmp(trLabelInput,'stim') && (all(ismember({'ci','av'},trCondsCell)) || ...
        all(ismember({'i','av'},trCondsCell)))
    error('mvpa:svmTrain:invalidInput',...
        ['Ambiguous training label for the specified training condition. ',...
        'Use a different label!']);
end

%% Checking...
% matlab version
if verLessThan('matlab', '8.3.0.532')
    error('mvpa:svmTrain:notSupportedMatlabVersion',...
        'MATLAB 2014a or newer is required! ');
end

% parallel pool 
currPool = gcp('nocreate');
if isempty(currPool)
    error('mvpa:svmTrain:missingParallelPool',...
        'Please open a parallel pool in order to run the function in parallel mode!');
else
    poolSize = currPool.NumWorkers;
end

% libsvm 
if isempty(regexp(path,'libsvm','once'))
    error('mvpa:svmTrain:missingToolbox','libsvm not found in path!');
end 

% parfor_progress
if progrMonitor
    if isempty(regexp(path,'parfor_progress','once'))
        error('mvpa:svmTrain:missingToolbox','parfor_progress not found in path!');
    end
end

%% Starting the timer
cStartFun = clock;
ds = datestr(now);
fprintf('\nsvmtrain%s%s\n%s\n',repmat(' ',1,72-length(ds)),ds,repmat('-',1,80));

%% Loading data and definig important variables
dataM = matfile(dirDataFile);
condDef = dataM.condDef;
misc = dataM.misc;
if any(strcmp(trMethod,{'sample-wise-beta','sample-wise-beta-rl','sample-wise-sm-beta'}))
    if strcmp(trLabelInput,'stim')
        infoStr = 'info_stim';
    elseif strcmp(trLabelInput,'resp')
        infoStr = 'info_resp';
    elseif strcmp(trLabelInput,'hand')
        error('mvpa:svmTrain:missingFunctionality',...
            'This analysis is not yet implemented!');
    end
else
    infoStr = 'info';
end
dataInfo = dataM.(infoStr);

%% Choosing the examples for training accroding to the defined conditions
if any(strcmp(trMethod,{'sample-wise-beta','sample-wise-beta-rl','sample-wise-sm-beta'})) && strcmp(trLabelInput,'resp')
    
    % In this case we don't really have stimulus conditions as all of them 
    % are represented in the glm, so we set the trCond manually and use all
    % the examples available. 
    if ~strcmp(trCond,'AV-ci-av')
        warning('mvpa:svmTrain:automaticOptionSetting',...
            'Changing training condition automatically to AV-ci-av!')
        trCond = 'AV-ci-av';
    end
    trIsExample = true(size(dataInfo,1),1);
    trLabel = 'resp';
    
else
    
    trIsExample = mvpa.selectexamples(trCond,condDef,dataInfo);
    if sum(trIsExample) == 0
        error('mvpa:svmTrain:requestedDataNotPresent',...
            'No examples of the specified condition were found in the dataset!');
    end
    trLabel = mvpa.selectlabel(trLabelInput,trCond);
    
end

%% Extracting the appropriate set of featureas and examples for training
% trExamples_feat is a 3D matrix with 
% size(1) = nFeatures
% size(2) = nTimePoints
% size(3) = nExamples
% Feature scaling is also done here within extractfeatures. 
trTimePointsSec = cellfun(@rdivide,trTimePoints,repmat({1000},size(trTimePoints)),'UniformOutput',false);
if strcmp(trMethod,'sample-wise-tf')
    trExamples_feat = mvpa.extractfeatures(dataM,trIsExample,trMethod,trLabelInput,trTimePointsSec,...
        'tfFreq',tfFreq,'tfType',tfType,'progrMonitor',progrMonitor);
elseif ismember(trMethod,{'sample-wise-bp','sample-wise-bp-avg'})
    trExamples_feat = mvpa.extractfeatures(dataM,trIsExample,trMethod,trLabelInput,trTimePointsSec,...
        'bpFreq',bpFreq,'progrMonitor',progrMonitor);
elseif ismember(trMethod,{'sample-wise-source','sample-wise-source-avg'})
    trExamples_feat = mvpa.extractfeatures(dataM,trIsExample,trMethod,trLabelInput,trTimePointsSec,...
        'trROI',trROI,'progrMonitor',progrMonitor);
else
    trExamples_feat = mvpa.extractfeatures(dataM,trIsExample,trMethod,trLabelInput,trTimePointsSec,...
        'progrMonitor',progrMonitor);
end
trExamples_info = dataInfo(trIsExample,:);

[nFeatures,nTimePoints,nExamples] = size(trExamples_feat);

%% Setting hyperparameters. 
switch svmType
    case 'cc'
        % Number of parameters required for the SVM
        nParam = 1;
        % For C-SVC, C is the cost parameter controlling the errors,
        % Here, log2(C) is given instead of C for computational convenience.
        % log2(C)
        if isempty(params), params{1} = 0; end
    case 'er'
        % Number of parameters required for the SVM
        nParam = 2;
        % For epsilon-SVR, C is the cost parameter controlling the errors,
        % e controls the width of the 'tube' which is fitted to the data.
        % Here, log2(C) and log2(e) is given for computational convenience.
        if isempty(params)
            params{1} = 0; % log2(C) or C = 1
            params{2} = -3.3219; % log2(e) or e = 0.1
        elseif numel(params) < 2
            error('mvpa:svmTrain:invalidInput',...
                'Either both hpyerparameters should be specified or none.');
        end
    case 'nc'
        % Number of parameters required for the SVM
        nParam = 1;
        % For nu-SVC, nu is the the upper and lower limit of errors and number
        % of support vectors respectively.
        % nu
        if isempty(params), params{1} = 0.5; end
    case 'nr'
        % Number of parameters required for the SVM
        nParam = 2;
        % For nuSVR, C is the cost parameter controlling the errors,
        % nu is the the upper and lower limit of errors and number of
        % support vectors respectively. Here, log2(C) is given instead of C 
        % for computational convenience.
        if isempty(params)
            params{1} = 0; % log2(C) or C = 1
            params{2} = 0.5; % nu
        elseif numel(params) < 2
            error('mvpa:svmTrain:invalidInput',...
                'Either both hpyerparameters should be specified or none! ');
        end
end

% Determining if gridsiearch is necessary
if nParam == 1
    if isscalar(params{1})
        doGridSearch = false;
    else
        doGridSearch = true;
        % Number of steps for the fine pass grid search.
        nStepsFine = 10;
    end
elseif nParam == 2
    if isscalar(params{1}) && isscalar(params{2})
        doGridSearch = false;
    else
        doGridSearch = true;
        % Number of steps for the fine pass grid search.
        nStepsFine = 10;
    end
end

%% Preparing variables for cross-validation
if strcmp(cvScheme,'kf')
    
    % Reseeding the random number generator and saving its status
    rng('shuffle');
    sRand = rng;
    
    % Preparing groupings for the cross-validation
    misc.k = k;
    misc.nKFrep = nKFrep;
    misc.nTimePoints = nTimePoints;
    groupings = mvpa.assigngroupings(cvScheme,trExamples_info,misc);
    % Cross-validation type specific input for the cross-validation 
    % function
    cvDetails.nKFrep = nKFrep;
    cvDetails.k = k;
    nModelsPerTimePoint = size(unique(groupings(:,1,1)),1)*size(groupings,2);
    % Choosing function according to the cross-validation scheme
    funCV = @mvpa.dokfoldcv;
    
elseif strcmp(cvScheme,'loso')
    
    % Preparing groupings for the cross-validation
    misc.nTimePoints = nTimePoints;
    groupings = mvpa.assigngroupings(cvScheme,trExamples_info,misc);
    % Cross-validation type specific input for the cross-validation 
    % function
    nModelsPerTimePoint = size(unique(groupings(:,1)),1);
    % Choosing function according to the cross-validation scheme
    funCV = @mvpa.dolosocv;
    
elseif strcmp(cvScheme,'loxo')
    error('mvpa:svmTrain:missingFunctionality',...
        'Cross-validation scheme ''loxo'' is not yet implemented! ');
end

% We need to generate the pseudoAV trials from A and V trials here
% and this will affect the number of examples in the dataset
if ~isempty(regexp(trCond,'pseudoAV','once'))
    % Saving the grouping corresponding the true A and V trials as
    % it will be needed for generalization
    groupings_AxorV = groupings;
    [trExamples_info,trExamples_feat,groupings] = ...
        mvpa.generatePseudoAvData(trExamples_info,trExamples_feat,groupings,trCond,condDef);
    % Updating the number of examples
    nExamples = size(trExamples_feat,3);
end

% The maximum possible number of support vectors over the folds is the
% number of examples in the largest of the training sets.
maxnSVs = nExamples - minnival(groupings(:,1,1));

% General input for the crossvalidation function
cvDetails.doGridSearch = doGridSearch;
cvDetails.nFeatures = nFeatures;
cvDetails.maxnSVs = maxnSVs;
cvDetails.nModelsPerTimePoint = nModelsPerTimePoint;
cvDetails.nParam = nParam;
cvDetails.params = params;
cvDetails.scMethod = scMethod;
cvDetails.svmType = svmType;


%% Training models
% Initializing progress monitor
cStartTr = clock;
fprintf('Training... \n');
if progrMonitor
    parfor_progress(nTimePoints);
end

% Assigning training labels
if ~ismember(trLabel,trExamples_info.Properties.VariableNames)
    error('mvpa:svmTrain:invalidVariableValue',...
        'The specified training label is not part of the dataset!');
else
    trLabs = trExamples_info.(trLabel);
end
% Checking training labels
if any(isnan(trLabs))
    error('mvpa:svmTrain:invalidVariableValue',...
        'Training labels can''t be NaNs!');
end

if any(ismember(svmType,{'cc','nc'}))
    nClassesMdl = max([2,numel(unique(trLabs))]);
else
    nClassesMdl = 2;
end

% Initializing arrays for data collection (padded with NaNs).
models = NaN(nFeatures+nClassesMdl+1,maxnSVs,nModelsPerTimePoint,nTimePoints);
scParam = NaN(nFeatures,2,nModelsPerTimePoint,nTimePoints);
acc = NaN(3,nModelsPerTimePoint,nTimePoints);
if doGridSearch
    [gridsCoarse,gridsFine,paramFine1,paramFine2,bestParCoarse,bestParFine] = deal(cell(nTimePoints,1));
    cvDetails.nStepsFine = nStepsFine;
end

% If the number of time pointsis greater or equal to the parallel pool 
% size, the parfor loop runs across the timepoints. 
if nTimePoints >= poolSize
    
    parfor iTimePoints = 1:nTimePoints
        % Extracting the actual time sample. This yields a matrix with number
        % of rows = number of examples, number of columns = number of features.
        actTimePointFeats = squeeze(trExamples_feat(:,iTimePoints,:))';
        actTimePointLabs = trLabs;
        actTimePointGroupings = [];
        if ismatrix(groupings)
            actTimePointGroupings = groupings(:,iTimePoints);
        elseif ndims(groupings) == 3
            actTimePointGroupings = groupings(:,:,iTimePoints);
        end
        % Collecting data.
        [models(:,:,:,iTimePoints),scParam(:,:,:,iTimePoints),acc(:,:,iTimePoints),cvMisc] = ...
            funCV(actTimePointFeats,actTimePointLabs,actTimePointGroupings,cvDetails);
        
        if doGridSearch
            gridsCoarse{iTimePoints} = cvMisc.gridsCoarse;
            gridsFine{iTimePoints} = cvMisc.gridsFine;
            paramFine1{iTimePoints} = cvMisc.paramFine1;
            if nParam == 2
                paramFine2{iTimePoints} = cvMisc.paramFine2;
            end
            bestParCoarse{iTimePoints} = cvMisc.bestParCoarse;
            bestParFine{iTimePoints} = cvMisc.bestParFine;
        end
        % Advancing Progress monitor if applicable
        if progrMonitor
            parfor_progress;
        end
    end
    
    
% If the number of time points is smaller than the parallel pool size, the 
% parfor loop runs across the parameter space of the gridsearch. 
else
    
    for iTimePoints = 1:nTimePoints
        
        % Extracting the actual time sample. This yields a matrix with number
        % of rows = number of examples, number of columns = number of features.
        actTimePointFeats = squeeze(trExamples_feat(:,iTimePoints,:))';
        actTimePointLabs = trLabs;
        if ismatrix(groupings)
            actTimePointGroupings = groupings(:,iTimePoints);
        elseif ndims(groupings) == 3
            actTimePointGroupings = groupings(:,:,iTimePoints);
        end
        % Collecting data.
        [models(:,:,:,iTimePoints),scParam(:,:,:,iTimePoints),acc(:,:,iTimePoints),cvMisc] = ...
            funCV(actTimePointFeats,actTimePointLabs,actTimePointGroupings,cvDetails);
        
        if doGridSearch
            gridsCoarse{iTimePoints} = cvMisc.gridsCoarse;
            gridsFine{iTimePoints} = cvMisc.gridsFine;
            paramFine1{iTimePoints} = cvMisc.paramFine1;
            if nParam == 2
                paramFine2{iTimePoints} = cvMisc.paramFine2;
            end
            bestParCoarse{iTimePoints} = cvMisc.bestParCoarse;
            bestParFine{iTimePoints} = cvMisc.bestParFine;
        end
        % Advancing Progress monitor if applicable
        if progrMonitor
            parfor_progress;
        end
        
    end
    
end

% Finalizing progress monitor.
if progrMonitor
    parfor_progress(0);
end
fprintf('Training elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(clock,cStartTr)/86400,'dd HH:MM:SS'));

%% Organizing and saving the output 
trainM = struct();
trainM.info.cv_scheme = cvScheme;
if strcmp(cvScheme,'kf')
    trainM.info.k = k;
    trainM.info.nKFrep = nKFrep;
end
trainM.info.params = params;
trainM.info.sc_method = scMethod;
trainM.info.subID = subID;
trainM.info.svm_type = svmType;
if strcmp(trMethod,'sample-wise-tf')
    trainM.info.tf_freq = tfFreq;
    trainM.info.tf_type = tfType;
elseif any(ismember(trMethod,{'sample-wise-bp','sample-wise-bp-avg'}))
    trainM.info.bp_freq = bpFreq;
elseif ismember(trMethod,{'sample-wise-source','sample-wise-source-avg'})
    trainM.info.tr_roi = trROI;
end
trainM.info.tr_cond = trCond;
trainM.info.tr_method = trMethod;
trainM.info.tr_timePoints = trTimePointsSec;
trainM.info.tr_elapsedTime = etime(clock,cStartFun);
trainM.info.tr_isExample = trIsExample;
if strcmp(cvScheme,'kf')
    trainM.info.tr_sRand = sRand;
end
[~,name,ext] = fileparts(dirDataFile);
trainM.info.tr_data_file = [name,ext];
trainM.info.tr_data_fs = misc.fs;
trainM.info.tr_label = trLabel;
trainM.info.sc_method = scMethod;

trainM.tr_examples = trExamples_info;

if doGridSearch
    trainM.tr_grids(1).pass = 'coarse';
    trainM.tr_grids(1).grids = cat(ndims(gridsCoarse{1})+1,gridsCoarse{:});
    trainM.tr_grids(1).param1 = params{1}';
    if nParam == 2
        trainM.tr_grids(1).param2 = params{2}';
    end
    trainM.tr_grids(1).bestPar = cat(ndims(bestParCoarse{1})+1,bestParCoarse{:});
    trainM.tr_grids(2).pass = 'fine';
    trainM.tr_grids(2).grids = cat(ndims(gridsFine{1})+1,gridsFine{:});
    trainM.tr_grids(2).param1 = cat(ndims(paramFine1{1})+1,paramFine1{:});
    if nParam == 2
        trainM.tr_grids(2).param2 = cat(ndims(paramFine2{1})+1,paramFine2{:});
    end
    trainM.tr_grids(2).bestPar = cat(ndims(bestParFine{1})+1,bestParFine{:});
end

trainM.tr_groupings = groupings;
if ~isempty(regexp(trCond,'pseudoAV','once'))
    trainM.tr_groupings_AxorV = groupings_AxorV;
end

trainM.tr_models = models;

trainM.tr_accuracies = acc;

if ~strcmp(scMethod,'Z-f')
    trainM.tr_scParam = scParam;
end

% Ordering the fields of the struct
trainM.info = orderfields(trainM.info);
trainM = orderfields(trainM); %#ok<*NASGU>

% Generating output file name.    
% Substring for the time samples
if size(trTimePoints,2) == 1
    if size(trTimePoints{1},2) > 1
        m = numel(trTimePoints{1});
        d = mean(diff(trTimePoints{1}));
        samplStr = sprintf('%d-m-%d-%d',mean(trTimePoints{1}),m,d);
    else
        samplStr = sprintf('%d',trTimePoints{1});
    end
elseif size(trTimePoints,2) > 1
    if any(cellfun(@numel,trTimePoints) > 1)
        m = numel(trTimePoints{1});
        d = mean(diff(trTimePoints{1}));
        trTimePoints = cellfun(@mean,trTimePoints);
        samplStr = sprintf('%d-%d-%d-m-%d-%d',min(trTimePoints),...
            round(mean(diff(trTimePoints))),max(trTimePoints),m,d);
    else
        trTimePoints = cell2mat(trTimePoints);
        samplStr = sprintf('%d-%d-%d',min(trTimePoints),round(mean(diff(trTimePoints))),...
            max(trTimePoints));
    end
end
% Substring for training label
 trLabelStr = trLabelInput;
 
% Special substring for the frequencies or ROI
if strcmp(trMethod,'sample-wise-tf')
    
   if ismember(tfType,{'f','ph'})
       s = cell(size(tfFreq{1}));
       for i = 1:size(tfFreq{1})
           if mod(tfFreq{1}(i),1), s{i} = '-%.1f'; else s{i} = '-%d'; end
       end
       strSlot = strrep(sprintf(['_%s' [s{:}]],tfType,tfFreq{1}),'.',',');
   else
       s = {cell(size(tfFreq{1})),cell(size(tfFreq{2}))};
       for i = 1:size(s,2)
           for j = 1:size(tfFreq{i},2)
               if mod(tfFreq{i}(j),1), s{i}{j} = '-%.1f'; else s{i}{j} = '-%d'; end
           end
       end
       strSlot = strrep(sprintf(['_f',[s{1}{:}],'_ph',[s{2}{:}]],tfFreq{1},tfFreq{2}),'.',',');      
   end

elseif any(ismember(trMethod,{'sample-wise-bp','sample-wise-bp-avg'}))
    strSlot = ['_bp-',bpFreq];
elseif ismember(trMethod,{'sample-wise-source','sample-wise-source-avg'})
    strSlot = ['_roi-',trROI];
else
    strSlot = '';
end

fileName = sprintf('%s_tr-%s_%s_%s%s_%s.mat',svmType,trCond,trLabelStr,...
    samplStr,strSlot,datestr(now,'yymmddHHMMSS'));
filePath = fullfile(dirAnalysis,fileName);

% Saving outputs. 
save(filePath,'-struct','trainM','-v7.3');
mvparesObj = mvpares(filePath);

end
