function genExamplesInfo = getInfoGenExamples(obj,varargin)
% Method for accessing the table of various information about the generalization examples
%
% USAGE:
%   genExamplesInfo = getInfoGenExamples(obj)
%   genExamplesInfo = getInfoGenExamples(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       cvFold (scalar/vector): the indices of cross-validation folds, 
%           default: all folds.
% OUTPUT:
%   genExamplesInfo (table): information about generalization examples

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
addRequired(p,'obj');
addParameter(p,'cvFold',[],@(x)validateattributes(x,...
    {'numeric'},{'vector','integer','positive'}));
parse(p,obj,varargin{:});
obj = p.Results.obj;
cvFold = p.Results.cvFold;

% Assign empty array and return if the result dataset is group level or it
% is not generalized. 
if strcmp(obj.level,'group')
    genExamplesInfo = [];
    warning('mvpares:getInfoGenExamples:datasetLevelMismatch',...
        ['The dataset''s level is ''group'', so it does ',...
        'not contain subject level information about generalization data']);
    return;
elseif strcmp(obj.state,'trained')
    genExamplesInfo = [];
    warning('mvpares:getInfoGenExamples:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end
% Checking cvFold
if isempty(cvFold)
    cvFold = 1:obj.getNcvFolds;
else
    if any(~ismember(cvFold,1:obj.getNcvFolds))
        error('mvpares:getInfoGenExamples:invalidInput',...
            'The input given for ''cvFold'' is invaid.');
    end
end

% Loading necessary data
cvScheme = obj.info.cv_scheme;
genExamples = obj.data.gen_examples;
if isprop(obj.data,'gen_groupings')
    groupings = obj.data.gen_groupings;
else
    groupings = obj.data.tr_groupings;
end

% Collecting true labels
if any(strcmp(obj.getTrGenCondMatch,{'full','partial','no'}))
    % If the generalization and training examples at least
    % partially overlap, we have to select the true labels
    % accroding to the cross-validation groupings.
    
    % Get rid of the redundant grouping info (it should be unifrom
    % across time)
    switch cvScheme
        case 'kf'
            % Check if the groupings are uniform across time (for
            % backward compatibility)
            sg = size(groupings);
            temp = mat2cell(groupings,ones(sg(1),1),ones(sg(2),1),sg(3));
            temp = cellfun(@squeeze,temp,'UniformOutput',false);
            temp = cellfun(@unique,temp,'UniformOutput',false);
            temp = cell2mat(cellfun(@numel,temp,'UniformOutput',false));
            if any(temp(:) > 1)
                error('mvpares:getInfoGenExamples:groupingInconsistent',...
                    'The CV grouping is not uniform across time points');
            end
            groupings = groupings(:,:,1);
        case 'loso'
            groupings = groupings(:,1);
    end
    
    % Finding the common group label if there is any (it is
    % a negative number, almost certainly -1)
    groupLabels = unique(groupings);
    commonGroupLabel = groupLabels(groupLabels < 0);
    if isempty(commonGroupLabel)
        commonGroupLabel = NaN;
    elseif ~isscalar(commonGroupLabel)
        error('mvpares:getInfoGenExamples:groupingInconsistent',...
            'The CV grouping contains more than one common label');
    end
    
    genExamplesInfo = table;
    switch cvScheme
        case 'kf'
            k = obj.info.k;
            for i = 1:numel(cvFold)
                iKFfold = mod(cvFold(i),k);
                if iKFfold == 0, iKFfold = k; end
                actGrouping = groupings(:,ceil(cvFold(i)/k));
                actFoldExamples = actGrouping == iKFfold | ...
                    actGrouping == commonGroupLabel;
                temp = genExamples(actFoldExamples,:);
                temp.fold = ones(size(temp,1),1)*cvFold(i);
                genExamplesInfo = cat(1,genExamplesInfo,temp);
            end
        case 'loso'
            for i = 1:numel(cvFold)
                actFoldExamples = groupings == cvFold(i) | ...
                    groupings == commonGroupLabel;
                temp = genExamples(actFoldExamples,:);
                temp.fold = ones(size(temp,1),1)*cvFold(i);
                genExamplesInfo = cat(1,genExamplesInfo,temp);
            end
    end
% elseif strcmp(obj.getTrGenCondMatch,'no')
%     % If there is no overlap between the training and
%     % generalization examples we use all the examples for
%     % each fold.
%     genExamplesInfo = repmat(genExamples,numel(cvFold),1);
%     temp = repmat(cvFold,size(genExamples,1),1);
%     genExamplesInfo.fold = temp(:);
else
    genExamplesInfo = [];
    warning('mvpares:getInfoGenExamples:undeterminedCase',...
        ['Could not determine the relationship between the training and ',...
        'generalization conditions.']);
    return;
end
s = size(genExamplesInfo);
% Moving variable 'fold' to the first place
genExamplesInfo = genExamplesInfo(:,[s(2),1:(s(2)-1)]);

end