function feat = extractfeatures(dataM,isExample,trMethod,trLabel,timePoints,varargin)
% Extracts the required set of features

%% Parsing input
p = inputParser;

% Valid values for certain inputs.
validBpFreq = {'d','t','a','b','gl','gh'};
validTrMethods = {'sample-wise','sample-wise-avg','sample-wise-beta',...
    'sample-wise-beta-rl','sample-wise-bp','sample-wise-bp-avg',...
    'sample-wise-rl','sample-wise-sc','sample-wise-mult','sample-wise-sm',...
    'sample-wise-sm-avg','sample-wise-sm-beta','sample-wise-sm-corr1', ...
    'sample-wise-sm-corr2','sample-wise-sm-avg-corr2','sample-wise-tf',...
    'sample-wise-source','sample-wise-source-avg'};

validTrLabels = {'disp','finger','hand','relV','resp','stim','task','sAhatInd',...
                 'sVhatInd','sHatComm','sAhatANDsVhat','LvsRa','LvsRv'};
validTfTypes = {'f','ph','fph'};
validROIs = {'A','IPS','V1-3'};
% Custom validation functions for certain inputs. 
    function checkTimePoints(x)
        if ~iscell(x) || ~isrow(x)
            error('mvpa:extractfeatures:invalidInput',...
                'tr_samples must be a cell array of one row!');
        elseif any(~cellfun(@isrow,x)) || numel(unique(cellfun(@numel,x))) > 1
            error('mvpa:extractfeatures:invalidInput',...
                'All elements of tr_timePoints must be row vectors of equal length!');
        end
    end

% Defining inputs.
addRequired(p,'dataM',@(x)validateattributes(x,{'matlab.io.MatFile'},{'nonempty'}));
addRequired(p,'isExample',@(x)validateattributes(x,{'logical'},{'column'}));
addRequired(p,'trMethod',@(x)any(validatestring(x,validTrMethods)));
addRequired(p,'trLabel',@(x)any(validatestring(x,validTrLabels)));
addRequired(p,'timePoints',@checkTimePoints);
addParameter(p,'bpFreq','',@(x)any(validatestring(x,validBpFreq)));
addParameter(p,'progrMonitor',true,@islogical);
addParameter(p,'tfFreq',{},@iscell);
addParameter(p,'tfType','',@(x)any(validatestring(x,validTfTypes)));
addParameter(p,'trROI','',@(x)any(validatestring(x,validROIs)));

% Parsing inputs.
parse(p,dataM,isExample,trMethod,trLabel,timePoints,varargin{:});

% Assigning inputs to variables.
dataM = p.Results.dataM;
isExample = p.Results.isExample;
trMethod = p.Results.trMethod;
trLabel = p.Results.trLabel;
timePoints = p.Results.timePoints;
bpFreq = p.Results.bpFreq;
progrMonitor = p.Results.progrMonitor;
tfFreq = p.Results.tfFreq;
tfType = p.Results.tfType;
trROI = p.Results.trROI;

% Custom checkings after parsing
if strcmp(trMethod,'sample-wise-tf') && any(cellfun(@isempty,{tfFreq,tfType}))
    error('mvpa:extractfeatures:missingInput',...
        'tfFreq and tfType must be specified!');
end

requiredVarsBP = {'tr_roi'};
if ismember(trMethod,{'sample-wise-source','sample-wise-source-avg'}) ...
        && any(ismember(requiredVarsBP,p.UsingDefaults))
    error('mvpa:extractfeatures:missingInput','''trROI'' must be specified!');
end

if strcmp(trMethod,'sample-wise-bp') && isempty(bpFreq)
    error('mvpa:extractfeatures:missingInput','bpFreq must be specified!');
end

%%
misc = dataM.misc;

if ismember(trMethod,{'sample-wise','sample-wise-avg','sample-wise-beta',...
        'sample-wise-beta-rl','sample-wise-sm-avg','sample-wise-sm-beta',...
        'sample-wise-sm-corr1','sample-wise-sm-corr2','sample-wise-sm-avg-corr2',...
        'sample-wise-rl','sample-wise-sc','sample-wise-mult','sample-wise-sm',...
        'sample-wise-source','sample-wise-source-avg'})
    
    if any(ismember(trMethod,{'sample-wise-beta','sample-wise-beta-rl','sample-wise-sm-beta'}))
        if strcmp(trLabel,'stim')
            featStr = 'feat_stim';
        elseif strcmp(trLabel,'resp')
            featStr = 'feat_resp';
        end
    elseif ismember(trMethod,{'sample-wise-source','sample-wise-source-avg'})
        featStr = sprintf('feat_%s',strrep(trROI,'-','_'));
    else
        featStr = 'feat';
    end
    dataFeatSize = size(dataM,featStr);
    nFeatures = dataFeatSize(1)*size(timePoints{1},2);
    nSamples = size(timePoints,2);
    nExamples = sum(isExample);
    feat = NaN(nFeatures,nSamples,nExamples);
    
    % Initializing progress monitor
    cStartExtr = clock;
    fprintf('Extracting features...\n');
    if progrMonitor
        parfor_progress(size(timePoints,2));
    end
    
    allFeats = dataM.(featStr);
    for i = 1:size(timePoints,2)
        
        actSample = mvpa.indsamplecustom(timePoints{i},dataFeatSize(2),misc.fs,misc.timeOnset);
        % Error if any specified training sample is out of the boundaries of the
        % data.
        if any(isnan(actSample))
            error('mvpa:extractfeatures:invalidSample',...
                'At least one specified sample is not part of the data set!');
        end
        
        % Select the actual time sample(s) and the appropriate examples
        temp = allFeats(:,actSample,isExample);
        
        % Reshaping the array to nFeatures x 1 x nExamples 
        s = size(temp);
        feat(:,i,:) = reshape(temp,s(1)*s(2),1,s(3));
        
        % Advancing progress bar
        if progrMonitor
            parfor_progress;
        end
        
    end
    
elseif ismember(trMethod,{'sample-wise-bp','sample-wise-bp-avg'})
    
    featStr = 'feat';
    dataFeatSize = size(dataM,featStr);
    nFeatures = dataFeatSize(1)*size(timePoints{1},2);
    nSamples = size(timePoints,2);
    nExamples = sum(isExample);
    feat = NaN(nFeatures,nSamples,nExamples);
    
    % Initializing progress monitor
    cStartExtr = clock;
    fprintf('Extracting features...\n');
    if progrMonitor
        parfor_progress(size(timePoints,2));
    end
    
    allFeats = dataM.(featStr);
    for i = 1:size(timePoints,2)
        
        actSample = mvpa.indsamplecustom(timePoints{i},dataFeatSize(3),misc.fs,misc.timeOnset);
        % Error if any specified training sample is out of the boundaries of the
        % data.
        if any(isnan(actSample))
            error('mvpa:extractfeatures:invalidSample',...
                'At least one specified sample is not part of the data set!');
        end
        
        % Error if any specified frequency is not part of the data
        if any(~ismember(bpFreq,misc.freqStr))
            error('mvpa:extractfeatures:invalidFrequencyBand',...
                ['At least one specified frequency band is not part of the ',...
                'data set!']);
        end
        
        actFreqInd = find(ismember(misc.freqStr,bpFreq));
        temp = allFeats(:,actFreqInd,actSample,isExample);
        
        % 1. Reshaping the array to nFeatures x 1 x nExamples
        s = size(temp);
        feat(:,i,:) = reshape(temp,s(1)*s(2)*s(3),1,s(4));
        
        % Advancing progress bar
        if progrMonitor
            parfor_progress;
        end
        
    end
    
elseif strcmp(trMethod,'sample-wise-tf')
    
    if ismember(tfType,{'f','ph'})
        
        if strcmp(tfType,'f')
            featStr = 'feat_f';
        else
            featStr = 'feat_ph';
        end
        dataFeatSize = size(dataM,featStr);
        nFeatures = dataFeatSize(1)*size(timePoints{1},2)*size(tfFreq{1},2);
        nSamples = size(timePoints,2);
        nExamples = sum(isExample);
        feat = NaN(nFeatures,nSamples,nExamples);
        
        % Initializing progress monitor
        cStartExtr = clock;
        fprintf('Extracting features...\n');
        if progrMonitor
            parfor_progress(size(timePoints,2));
        end
        
        allFeats = dataM.(featStr);
        for i = 1:size(timePoints,2)
            
            actSample = mvpa.indsamplecustom(timePoints{i},dataFeatSize(3),misc.fs,misc.timeOnset);
            % Error if any specified training sample is out of the boundaries of the
            % data.
            if any(isnan(actSample))
                error('mvpa:extractfeatures:invalidSample',...
                    'At least one specified sample is not part of the data set!');
            end
            
            actFreq = tfFreq{1};
            % Error if any specified frequency is not part of the data
            if any(~ismember(actFreq,misc.frequencies))
                error('mvpa:extractfeatures:invalidFrequency',...
                    'At least one specified frequency is not part of the data set!');
            end
            actFreqInd = find(ismember(misc.frequencies,actFreq));
            temp = allFeats(:,actFreqInd,actSample,isExample); %#ok<*FNDSB>
                        
            % 1. Reshaping the array to nFeatures x 1 x nExamples
            s = size(temp);
            feat(:,i,:) = reshape(temp,s(1)*s(2)*s(3),1,s(4));
            
            % Advancing progress bar
            if progrMonitor
                parfor_progress;
            end
            
        end
        
    elseif strcmp(tfType,'fph')
        
        featStr = {'feat_f','feat_ph'};
        dataFeatSize = size(dataM,featStr{1});
        nFeatures = 0;
        for i = 1:size(featStr,2)
            inc = dataFeatSize(1)*size(timePoints{1},2)*size(tfFreq{i},2);
            nFeatures = nFeatures + inc;
        end
        nSamples = size(timePoints,2);
        nExamples = sum(isExample);
        feat = NaN(nFeatures,nSamples,nExamples);
        
        % Initializing progress monitor
        cStartExtr = clock;
        fprintf('Extracting features...\n');
        if progrMonitor
            parfor_progress(size(timePoints,2));
        end
        % This needs further optimization
        parfor i = 1:size(timePoints,2)
            
            actSample = mvpa.indsamplecustom(timePoints{i},dataFeatSize(3),misc.fs,misc.timeOnset);
            % Error if any specified training sample is out of the boundaries of the
            % data.
            if any(isnan(actSample))
                error('mvpa:extractfeatures:invalidSample',...
                    'At least one specified sample is not part of the data set!');
            end
            
            temp = cell(1,size(featStr,2));
            for j = 1:size(featStr,2)
                
                actFreq = tfFreq{j};
                % Error if any specified frequency is not part of the data
                if any(~ismember(actFreq,misc.frequencies))
                    error('mvpa:extractfeatures:invalidFrequency',...
                        'At least one specified frequency is not part of the data set!');
                end
                actFreqInd = find(ismember(misc.frequencies,actFreq));
                temp{j} = dataM.(featStr{j})(:,actFreqInd,actSample,:);
                temp{j} = temp{j}(:,:,:,isExample);
                
                % 1. Reshaping the array to nFeatures x 1 x nExamples
                s = size(temp{j});
                temp{j} = reshape(temp{j},s(1)*s(2)*s(3),1,s(4));
                
            end
            % The two types of features - power and phase - are scaled
            % separately and now concatenated to one featre set. 
            feat(:,i,:) = cat(1,temp{:});
            
            % Advancing progress bar
            if progrMonitor
                parfor_progress;
            end
        end
    
    end
    
end

% Finalizing progress monitor 
if progrMonitor
    parfor_progress(0);
end
fprintf('Feature extraction elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(clock,cStartExtr)/86400,'dd HH:MM:SS'));

end

