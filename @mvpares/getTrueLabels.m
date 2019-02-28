function trueLabels = getTrueLabels(obj,varargin)
% Method for getting true labels
%
% USAGE:
%   trueLabels = getTrueLabels(obj)
%   trueLabels = getTrueLabels(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       cvFold (scalar/vector): the indices of cross-validation folds. 
%           default: all folds
% OUTPUT:
%   trueLabels (array): predicted labels, size depends on
%       the particular mvpa result dataset, and input settings.

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

% Assign empty array and return if the dataset is of group level 
if strcmp(obj.level,'group')
    trueLabels = [];
    warning('mvpares:getTrueLabels:datasetLevelMismatch',...
        ['The dataset''s level is ''group'', so it does ',...
        'not contain predicted labels.']);
    return;
end
% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(obj.state,'trained')
    trueLabels = [];
    warning('mvpares:getTrueLabels:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end
% Check if there is a generalization label
if isempty(obj.info.gen_label)
    trueLabels = [];
    warning('mvpares:getTrueLabels:requestedDataNotPresent',...
        'The dataset does not specifies a generalization label.');
    return;
end
% Checking cvFold
if isempty(cvFold)
    cvFold = 1:obj.getNcvFolds;
else
    if any(~ismember(cvFold,1:obj.getNcvFolds))
        error('mvpares:getTrueLabels:invalidInput',...
            'The input given for ''cvFold'' is invaid.');
    end
end

labelToUse = obj.info.gen_label;

temp = obj.getInfoGenExamples('cvFold',cvFold);

% Checking if the generalization label is part of the
% generailzation examples table. Merging datasets with different
% lables do not have one generalization label.
if ~ismember(labelToUse,temp.Properties.VariableNames)
    trueLabels = [];
    warning('mvpares:getTrueLabels:requestedDataNotPresent',...
        ['The generalization examples do not contain the variable ' ...
         '%s. Returning'],labelToUse);
    return;
end

trueLabels = table2array(temp(:,labelToUse));

end