function genTimePoints = getGenTimePoints(obj)
% Method for getting the generalization time points
% 
% USAGE:
%   genTimePoints = getGenTimePoints(obj)
% INPUT:
%   obj (object): mvpares object
% OUTPUT:
%   genTimePoints (vector/matrix): generalization time points. If the
%       dataset is generalized to training time it is a vector. If the
%       dataset is generalized to training time x training time it is a
%       matrix.
%
% Copyright(C) 2016, Mate Aller

% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(obj.state,'trained')
    genTimePoints = [];
    warning('mvpares:getGenTimePoints:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end

genTimePoints = obj.getTrTimePoints;
if strcmp(obj.info.gen_time,'tr_x_tr')
    genTimePoints = repmat(genTimePoints',size(genTimePoints,1),1);
end

end