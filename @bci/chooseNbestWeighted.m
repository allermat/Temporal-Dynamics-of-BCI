function idx = chooseNbestWeighted(param,perf,n,thres)
% Find the N best performing parameter combinations weighted by their Eucledian distance
% 
% DETAILS: First find the thres number of best performing parameter
%   combinations. Then order these by their Eucledian distance from the 
%   best perfoming and sample n from them so the samples are at equal
%   distance from each other and cover the whole selected range. 
% USAGE: 
%   idx = chooseNbestWeighted(param,perf,n,thres)
% INPUT:
%   param (numeric matrix)
%   perf (numeric vector)
%   n (numeric scalar)
%   thres (numeric scalar)
% OUTPUT: 
%   idx (numeric vecor): indices of the n choosen parameter combintations

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;

addRequired(p,'param',@(x) validateattributes(x,{'numeric'},{'nonempty'}));
addRequired(p,'perf',@(x) validateattributes(x,{'numeric'},{'vector',...
    'numel',size(param,1)}));
addRequired(p,'n',@(x) validateattributes(x,{'numeric'},{'scalar',...
    '<=',size(param,1)}));
addRequired(p,'thres',@(x) validateattributes(x,{'numeric'},{'scalar',...
    '<=',size(param,1)}));
parse(p,param,perf,n,thres);
param = p.Results.param;
perf = p.Results.perf;
n = p.Results.n;
thres = p.Results.thres;

% Sort the performances of parameter combinations
[~,idx] = sort(perf);
% Choose the best parameter combination(i.e. lowest logLike)
% and the following thres number of candidate combinations
bestParam = param(idx(1),:);
candidParam = param(idx(1:thres),:);
idxCandid = idx(1:thres);
% Compute the Eucledian distance between each candidate parameter vector
% and the best parameter combination
diff = candidParam - repmat(bestParam,size(candidParam,1),1);
dist = cellfun(@norm,mat2cell(diff,ones(size(diff,1),1),size(diff,2)));
% Sorting distances to ascending order
[~,distSortIdx] = sort(dist);
% Choosing n parameter vectors from the candidates sampling
idx = idxCandid(distSortIdx(round(linspace(1,thres,n))));


end

