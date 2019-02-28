function maxNum = maxnival(v)
% Retruns the maximum number of identical values in a vector. 
% 
% DETAILS:
%   Infs, -Infs and NaNs are automatically ignored. 
% 
% USAGE: 
%   maxNum = maxnival(v)
%
% INPUT:
%   v: vector of numbers. 
% 
% OUTPUT: 
%   maxNum: the maximum number of identical values in v. 

%% Parse input

p = inputParser;

addRequired(p,'v',@isvector);

parse(p,v);

v = p.Results.v;

%%
% Converting v to row vector if it was colum. 
if iscolumn(v)
    v = v';
end

% Removing Infs, -Infs and NaNs
v = v(isfinite(v));
u = unique(v); 

maxNum = 1;

if numel(u) == numel(v)
    return;
end

for i = 1:numel(u)
    
    num = sum(v == u(i));
    
    if num > maxNum
        maxNum = num;
    end
    
end


end

