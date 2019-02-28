function minNum = minnival(v)
% Retruns the minimum number of identical values in a vector. 
% 
% DETAILS:
%   Infs, -Infs and NaNs are automatically ignored. 
% 
% USAGE: 
%   minNum = minnival(v)
%
% INPUT:
%   v: vector of numbers. 
% 
% OUTPUT: 
%   minNum: the minimum number of identical values in v. 
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
       
if numel(u) == numel(v)
    minNum = 1;
    return;
elseif numel(u) == 1
    minNum = numel(v);
    return;
end

minNum = numel(v);
for i = 1:numel(u)
    
    num = sum(v == u(i));
    
    if num < minNum
        minNum = num;
    end
    
end

end

