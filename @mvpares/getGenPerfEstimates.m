function estimates = getGenPerfEstimates(obj,varargin)
% Method for getting prediction performance estimates
%
% USAGE:
%   estimates = getGenPerfEstimates(obj)
%   estimates = getGenPerfEstimates(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
% OUTPUT:
%   estimates (struct): prediction performance estimates as fields of a 
%       structure array. Each field contains a time series of estimates, 
%       or a time x time matrix of estimates, depending on the input
%       settings.
% 
% Copyright(C) 2016, Mate Aller

% Parsing input
p = inputParser;
validGenTimes = {'tr','tr_x_tr'};
addRequired(p,'obj');
addParameter(p,'genTime','tr',@(x)any(validatestring(x,validGenTimes)));
parse(p,obj,varargin{:});
obj = p.Results.obj;
genTime = p.Results.genTime;

% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(obj.state,'trained')
    estimates = [];
    warning('mvpares:getGenPerfEstimates:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end

% Assign empty array and return if the result dataset does not contain
% prediction performance estimates
if ~ismember('gen_perfEstimates',obj.who)
    estimates = [];
    warning('mvpares:getGenPerfEstimates:reqestedDataNotPresent',...
        'The dataset does not contain prediction performance estimates. ');
    return;
end

if strcmp(genTime,'tr_x_tr')
    % Check if genTime is valid given the result dataset
    if isvector(obj.getGenTimePoints)
        estimates = [];
        warning('mvpares:getGenPerfEstimates:reqestedDataNotPresent',...
            'The dataset is not generalized time x time! ');
        return;
    else
        estimates = obj.data.gen_perfEstimates;
    end
elseif strcmp(genTime,'tr')
    estimates = obj.data.gen_perfEstimates;
    if ~isvector(obj.getGenTimePoints)
        estimates = structfun(@diagcustom,estimates,'UniformOutput',false);
    end
end

end

function out = diagcustom(a)
% Custom function for taking the diagonal of matrices of 3D arrays
% For matrices it is the regular diag function.
% For 3D arrays the diag function is evaluated separately for each layer
% (third dimension). The output is a matrix, each column corresponding to
% each layer. 

if ismatrix(a)
    out = diag(a);
elseif ndims(a) == 3
    s = size(a);
    out = NaN(s(1),s(3));
    for i = 1:s(3)
        out(:,i) = diag(a(:,:,i));
    end
else
    error('mvpares:getGenPerfEstimates:diagcustom:invalidInput',...
        'Input must be matrix, or 3D array.');
end

end
