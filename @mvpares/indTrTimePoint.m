function idx = indTrTimePoint(obj,timePoint)
% Method for finding the index of a particular training time point
% 
% USAGE:
%   idx = indTrTimePoint(obj,timePoint)
% INPUT:
%   obj (object): mvpares object
%   timePoint (scalar): time point of interest in seconds
% OUTPUT:
%   idx (scalar): index of the time point of interest, [] if not found. 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

p = inputParser;

addRequired(p,'obj');
addRequired(p,'timePoint',@(x) validateattributes(x,{'numeric'},{'scalar'}));
parse(p,obj,timePoint);
obj = p.Results.obj;
timePoint = p.Results.timePoint;

time = obj.getTrTimePoints;
% eps(0.5) is the tolerance level for floating point comparison
idx = find(abs(time-timePoint) <= eps('double'));

end

