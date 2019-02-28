function mvparesObj = addAVmodelCorrelations(mvparesObj,varargin)
% Method for computing correleations of AV weights
% 
% USAGE:
%   mvparesObj = addAVmodelCorrelations(mvparesObj)
%   mvparesObj = addAVmodelCorrelations(mvparesObj,'Name',Value)
% INPUT:
%   Required:
%       mvparesObj (object): mvpares object
%   'Name'-Value arguments:
%       pathFmriFile (string): full path to the file containing fMRI data
%           for correlation
%       pathBehavFile (string): full path to the file containing 
%           behavioural data for correlation
% OUTPUT:
%   mvparesObj (object): mvpares object with AV model correlations

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
wrapexist = @(x) exist(x,'file');
addRequired(p,'mvparesObj',@(x) isa(x,'mvpares') && x.isvalid);
addParameter(p,'pathFmriFile','',wrapexist);
addParameter(p,'pathBehavFile','',wrapexist);

parse(p,mvparesObj,varargin{:});

mvparesObj = p.Results.mvparesObj;
pathFmriFile = p.Results.pathFmriFile;
pathBehavFile = p.Results.pathBehavFile;

% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(mvparesObj.state,'trained')
    mvparesObj = [];
    warning('mvpa:addAVmodelCorrelations:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end

% Starting the timer and printing details
cStart = clock;
ds = datestr(now);
printFnTitle(80,'addAVmodelCorrelations',ds)
fprintf('Estimating... \n');
% Selecting conditions
condsOfInterest = {'R_D_a','R_D_v','R_d_a','R_d_v',...
    'r_D_a','r_D_v','r_d_a','r_d_v'};
nTrTimePoints = mvparesObj.getSizeTrTime;
nTrTimePoints = nTrTimePoints(1);

for i = 1:2
    % Struct for collecting correlation data
    correlations = struct();
    % Loading AV model weights
    if i == 1
        genTime = mvparesObj.getInfo.gen_time;
        weights = mvparesObj.getAVmodelWeights('genTime',genTime);
    else
        weights = mvparesObj.getAVmodelWeights('genTime','tr','smooth',true);
    end
    % Preparing data for correlation
    fieldNames = fieldnames(weights);
    reqFieldIdx = ismember(fieldNames,strcat(condsOfInterest,'_wav'));
    temp = struct2cell(weights);
    temp = temp(reqFieldIdx);
    dimToCat = sum(size(temp{1}) > 1)+1;
    temp = cat(dimToCat,temp{:});
    if dimToCat == 2
        % This results a cell array of nTrTimePoints x 1
        toCorr = mat2cell(temp,ones(nTrTimePoints,1),numel(condsOfInterest));
        toCorr = cellfun(@squeeze,toCorr,'UniformOutput',false);
    elseif dimToCat == 3
        % This results a cell array of nTrTimePoints x nTrTimePoints
        toCorr = mat2cell(temp,ones(nTrTimePoints,1),ones(nTrTimePoints,1),numel(condsOfInterest));
        toCorr = cellfun(@squeeze,toCorr,'UniformOutput',false);
    else
        mvparesObj = [];
        warning('mvpa:addAVmodelCorrelations:unsupportedDimension',...
            ['The number of dimensions in the dataset is not supported. ',...
            'Please check dataset. Returning. ']);
        return;
    end
    
    % Correlating across time
    if size(toCorr,2) == 1
        toCorrTxT = repmat(toCorr,1,nTrTimePoints);
    else
        toCorrTxT = repmat(toCorr(logical(eye(size(toCorr,1)))),1,nTrTimePoints);
    end
    correlations.cc_acrossTime = cellfun(@circ_corrcc,toCorrTxT,toCorrTxT');
    
    % Correlating behavioural results
    if isempty(pathBehavFile)
        warning('mvpa:addAVmodelCorrelations:missingInput',...
            'No behavioural dataset was specified, skipping correlating with behavioural.');
    else
        temp = load(pathBehavFile);
        weights_behav = temp.wav;
        
        % Checking if all required conditions are present
        if any(~ismember(strcat(condsOfInterest,'_wav'),weights_behav.Properties.VariableNames))
            warning('mvpa:addAVmodelCorrelations:conditionMismatch',...
                ['At least one of the requested conditions is not found in the ',...
                'behavioural condition list. Skipping correlating with behavioural.']);
        end
        % Extracting
        temp = weights_behav(weights_behav.groupings == 'all',...
            strcat(condsOfInterest,'_wav'));
        temp = table2array(temp)';
        toCorrBehav = repmat({temp},size(toCorr));
        correlations.cc_behav = cellfun(@circ_corrcc,toCorr,toCorrBehav);
    end
    
    % Correlating fMRI results
    if isempty(pathFmriFile)
        warning('mvpa:addAVmodelCorrelations:missingInput',...
            'No fMRI dataset was specified, skipping correlating with fMRI.');
    else
        temp = load(pathFmriFile);
        tempFields = fieldnames(temp);
        weights_fmri = temp.(tempFields{1});
        % First sort the rows according to {'ROI','relV','discrepancy','task'},
        % then convert condition labels and select examples.
        weights_fmri = sortrows(weights_fmri,{'ROI','relV','discrepancy','task'});
        rois = categories(weights_fmri.ROI);
        temp = weights_fmri(weights_fmri.ROI == rois{1},{'relV','discrepancy','task'});
        temp = table2cell(temp);
        temp = cellfun(@char,temp,'UniformOutput',false);
        conds = strcat(temp(:,1),'_',temp(:,2),'_',temp(:,3));
        % Checking if all required conditions are present
        if any(~ismember(condsOfInterest,conds))
            warning('mvpa:addAVmodelCorrelations:conditionMismatch',...
                ['At least one of the requested conditions is not found in the ',...
                'fMRI condition list. Skipping correlating with fMRI.']);
        end
        [~,condIdx] = ismember(condsOfInterest,conds);
        
        for iRoi = 1:numel(rois)
            temp = weights_fmri.wav(weights_fmri.ROI == rois{iRoi});
            temp = temp(condIdx);
            toCorrFmri = repmat({temp},size(toCorr));
            fieldName = strcat('cc_fmri_',rois{iRoi});
            correlations.(fieldName) = cellfun(@circ_corrcc,toCorr,toCorrFmri);
        end
        
    end
    
    % Setting the dataset object to writable if it is not
    if ~mvparesObj.writable, mvparesObj.setWritable(true); end
    if i == 1
        fieldName = 'gen_AVmodelCorrelations';
    else
        fieldName = 'gen_AVmodelCorrelations_s';
    end
    if ismember(fieldName,fieldnames(mvparesObj.data))
        warning('mvpa:addAVmodelCorrelations:overwriteField',...
            ['The field %s already exists in the mvpa result dataset, '...
            'it will be overwritten.'],fieldName);
    end
    mvparesObj.data.(fieldName) = orderfields(correlations);
    mvparesObj.setWritable(false);
    
end

% Finishing timer and printing elapsed time
fprintf('Estimation elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

end

