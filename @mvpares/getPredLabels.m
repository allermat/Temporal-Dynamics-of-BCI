function predLabels = getPredLabels(obj,varargin)
% Method for getting predicted labels
%
% USAGE:
%   predLabels = getPredLabels(obj)
%   predLabels = getPredLabels(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       trTimePoints (scalar/vector): training time points
%           (default: all available time points)
%       genTime (string): indicating the generalization time. Possible
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'
%       cvFold (scalar/vector): the indices of cross-validation folds. 
%           default: all folds
% OUTPUT:
%   predLabels (array): predicted labels, size depends on
%       the particular mvpa result dataset, and input settings.

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validGenTimes = {'tr','tr_x_tr'};
addRequired(p,'obj');
addParameter(p,'cvFold',[],@(x)validateattributes(x,...
    {'numeric'},{'vector','integer','positive'}));
addParameter(p,'trTimePoints',obj.getTrTimePoints,@(x)validateattributes(x,...
    {'numeric'},{'vector','finite'}));
addParameter(p,'genTime','tr',...
    @(x)any(validatestring(x,validGenTimes)));
parse(p,obj,varargin{:});
obj = p.Results.obj;
cvFold = p.Results.cvFold;
trTimePoints = p.Results.trTimePoints;
genTime = p.Results.genTime;

% Assign empty array and return if the dataset is of group level 
if strcmp(obj.level,'group')
    predLabels = [];
    warning('mvpares:getPredLabels:datasetLevelMismatch',...
        ['The dataset''s level is ''group'', so it does ',...
        'not contain predicted labels.']);
    return;
end
% ...or it is not generalized. 
if strcmp(obj.state,'trained')
    predLabels = [];
    warning('mvpares:getPredLabels:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end
% Check if genTime is valid given the result dataset
if isvector(obj.getGenTimePoints) && strcmp(genTime,'tr_x_tr')
    predLabels = [];
    warning('mvpares:getPredLabels:reqestedDataNotPresent',...
        'The dataset is not generalized time x time! ');
    return;
end
% Checking requested training time points
if any(~ismember(trTimePoints,obj.getTrTimePoints))
    error('mvpares:getPredLabels:invalidInput',...
        ['The dataset does not contain at least one of the ',...
        'specfied training time points! ']);
end
% Checking cvFold
if isempty(cvFold)
    cvFold = 1:obj.getNcvFolds;
else
    if any(~ismember(cvFold,1:obj.getNcvFolds))
        error('mvpares:getPredLabels:invalidInput',...
            'The input given for ''cvFold'' is invaid.');
    end
end

sizePredLab = size(obj.data,'gen_predlabels');
% Indices for training and generalization time points
trTimeIdx = find(ismember(obj.getTrTimePoints,trTimePoints));
% Number of generalization time points per training time point
if strcmp(genTime,'tr')
    nGenTimePointsPerTrTimePoint = 1;
else
    nGenTimePointsPerTrTimePoint = sizePredLab(4);
end
% Loading data in one piece
predLabelsAll = obj.data.gen_predlabels;
% Pre-allocating output array
nExamplesOverall = size(obj.getInfoGenExamples('cvFold',cvFold),1);
predLabels = NaN(nExamplesOverall,numel(trTimeIdx),nGenTimePointsPerTrTimePoint);
startIdx = 1;
for i = 1:numel(cvFold)
    % Accessing data from the respective fold
    if numel(sizePredLab) < 4
        temp = squeeze(predLabelsAll(:,cvFold(i),trTimeIdx));
    else
        temp = squeeze(predLabelsAll(:,cvFold(i),trTimeIdx,:));
        % Selecting just the training time points
        if strcmp(genTime,'tr')
            idx = logical(eye(sizePredLab(3)));
            idx = idx(trTimeIdx,:);
            idx = shiftdim(idx,-1);
            idx = repmat(idx,sizePredLab(1),1,1);
            st = num2cell(size(temp));
            temp = temp(idx);
            temp = reshape(temp,st{1:2});
        end
    end
    % Getting rid of padded NaNs
    nExamplesActFold = size(obj.getInfoGenExamples('cvFold',cvFold(i)),1);
    temp = temp(1:nExamplesActFold,:,:,:);
    % Filling up the output array
    endIdx = startIdx+size(temp,1)-1;
    predLabels(startIdx:endIdx,:,:) = temp;
    startIdx = endIdx+1;
end

end