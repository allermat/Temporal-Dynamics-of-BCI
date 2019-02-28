function setPlotTimeWin(obj,timeWin)
% Method for setting the time window for plotting
%
% USAGE:
%   setPlotTimeWin(obj,timeWin)
% INPUT:
%   obj (object): mvpares object
%   timeWin (vector): time window for plotting in seconds
% 
% Copyright(C) 2016, Mate Aller

% Parsing input
p = inputParser;
addRequired(p,'obj');
addRequired(p,'timeWin',@(x) validateattributes(x,{'numeric'},...
    {'vector','numel',2,'increasing'}));
parse(p,obj,timeWin);
obj = p.Results.obj;
timeWin = p.Results.timeWin;

% Checking if specified time window is valid
if any(~ismember(timeWin,obj.getTrTimePoints))
    error('mvpares:setPlotTimeWin:invalidInput',...
        ['The dataset does not contain at least one of the specfied time ',... 
        'points! ']);
end

obj.plotTimeWin = timeWin;

end

