function out = getTrGenCondMatch(obj)
% Method for determining the relationship between the training and generalization conditions.
% 
% USAGE:
%   out = getTrGenCondMatch(obj)
% INPUT:
%   obj (object): mvpares object
% OUTPUT:
%   out (string): the relationship between the training and generalization 
%       conditions. 
%       Possible values: 'full' - complete match, 'partial' - partial match
%       'no' - disjoint sets of training and generalization examples

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

if strcmp(obj.level,'group')
    out = [];
    warning('mvpares:getTrGenCondMatch:datasetLevelMismatch',...
        ['The dataset''s level is ''group'', so it does ',...
        'not contain predicted labels.']);
    return;
end

if strcmp(obj.state,'trained')
    out = [];
    warning('mvpares:getTrGenCondMatch:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end
trIsExample = obj.info.tr_isExample;
genIsExample = obj.info.gen_isExample;
% This check is just for backward compatibility reasons
if ismember(obj.who,'tr_data_file')
    trDataFile = obj.info.tr_data_file;
    genDataFile = obj.info.gen_data_file;
else
    % Comparing two empty strings should give true
    trDataFile = '';
    genDataFile = '';
end

if ~strcmp(trDataFile,genDataFile)
    % If the training and generalization data files are different
    out = 'no';
elseif all(~ismember(find(genIsExample),find(trIsExample)))
    % If the training and generalization example sets are disjoint
    out = 'no';
elseif all(ismember(find(genIsExample),find(trIsExample)))
    % If the generalization example set matches the training example set
    % completely
    out = 'full';
else
    % If the generalization example set contains examples which are not
    % part of the training example set
    out = 'partial';
end

end