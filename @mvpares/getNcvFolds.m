function out = getNcvFolds(obj)
% Method for getting the number of cross-validation folds

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

if strcmp(obj.level,'group')
    warning('mvpares:getNcvFolds:datasetLevelMismatch',...
        ['The dataset''s level is ''group'', so it does not contain',...
        'this piece of information.']);
    out = [];
    return;
else
    if strcmp(obj.state,'trained')
        s = size(obj.data,'tr_models');
        out = s(3);
    elseif strcmp(obj.state,'trained_and_generalized')
        s = size(obj.data,'gen_predlabels');
        out = s(2);
    end
end
end