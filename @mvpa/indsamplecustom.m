function res = indsamplecustom(t,nsample,fs,starttime)
% Returns the sample closest to some time point in the specified time point set. 
% 
% USAGE: res = indsamplecustom(t,nsample,fs,starttime)
% 
% INPUT: 
%   t: row vector of time points in seconds. 
%   nsample: number of time samples in the time point set. 
%   fs: sampling frequency in Hz.
%   starttime: time of the first time sample in s. 
%
% OUTPUT:
%   res: row vector of indicies corresponding to the specifed time points.
%

%% Parsing input
p = inputParser;

% Defining inputs.
addRequired(p,'t',@(x)validateattributes(x,{'numeric'},{'row','finite','nonnan'}));
addRequired(p,'nsample',@(x)validateattributes(x,{'numeric'},{'scalar','finite','nonnan','positive'}));
addRequired(p,'fs',@(x)validateattributes(x,{'numeric'},{'scalar','finite','nonnan','positive'}));
addRequired(p,'starttime',@(x)validateattributes(x,{'numeric'},{'scalar','finite','nonnan'}));

% Parsing inputs.
parse(p,t,nsample,fs,starttime);

% Assigning inputs to variables.
t = p.Results.t;
nsample = p.Results.nsample;
fs = p.Results.fs;
starttime = p.Results.starttime;

%%
res = NaN(1,length(t));
T = (((1:nsample)-1)/fs)+starttime;
for i = 1:length(t)   
    [m,res(i)] = min(abs(T-t(i)));
    if m > (1/fs)
        warning('Could not find an index matching the requested time %d sec', t(i));
        res(i) = NaN;
    end
end

end

