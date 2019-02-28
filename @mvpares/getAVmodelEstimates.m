function estimates = getAVmodelEstimates(obj,varargin)
% Method for accessing AV model estimates
%
% USAGE:
%   estimates = getAVmodelEstimates(obj)
%   estimates = getAVmodelEstimates(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       smooth (logical): whether to smooth the time course of AV estimates
% OUTPUT:
%   estimates (struct): estimates as fields of the array. Each field
%       contains a time series of estimates of a certain condition, or a
%       time x time matrix of estimates depending on input settings. 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validGenTimes = {'tr','tr_x_tr'};
addRequired(p,'obj');
addParameter(p,'genTime','tr',@(x)any(validatestring(x,validGenTimes)));
addParameter(p,'smooth',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
parse(p,obj,varargin{:});
obj = p.Results.obj;
genTime = p.Results.genTime;
smooth = p.Results.smooth;

% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(obj.state,'trained')
    estimates = [];
    warning('mvpares:getAVmodelEstimates:datasetStateMismatch',...
        ['The dataset''s state is ''traned'', so it does ',...
        'not contain generalization data']);
    return;
end

% Assign empty array and return if the result dataset does not contain
% AV model estimates
if ~ismember('gen_AVmodelEstimates',obj.who)
    estimates = [];
    warning('mvpares:getAVmodelEstimates:reqestedDataNotPresent',...
        'The dataset does not contain AV model estimates. ');
    return;
end

% Loading estimates
estimates = obj.data.gen_AVmodelEstimates;

if strcmp(genTime,'tr_x_tr')
    % Check if genTime is valid given the result dataset
    if isvector(obj.getGenTimePoints)
        estimates = [];
        warning('mvpares:getAVmodelEstimates:reqestedDataNotPresent',...
            'The dataset is not generalized time x time! ');
        return;
    elseif smooth
        estimates = [];
        warning('mvpares:getAVmodelEstimates:invalidSettings',...
            ['Can''t smooth if time x time generalized estimates ',...
            'are required! ']);
        return;
    end
elseif strcmp(genTime,'tr')
    % Extracting the diagonal if just the training timepoints are needed
    % and the data is time x time generalized
    if ~isvector(obj.getGenTimePoints)
        estimates = structfun(@diag,estimates,'UniformOutput',false);
    end
    
    if smooth
        if strcmp(obj.level,'group')
            % For group level AV weights we save the smoothed time courses
            % as the smoothing happens on the subject level
            estimates = obj.data.gen_AVmodelEstimates_s;
        else
            % Smoothing the time courses moving average window of 20 ms
            fs = obj.getFsample;
            mavgWinSec = 0.02;
            mavgWinSample = round(mavgWinSec*fs);
            % Kernel for smoothing trials with the given time wintow
            kern = ones(mavgWinSample,1)./mavgWinSample;
            % Moving average function
            funMavg = @(x) conv(x,kern,'same');
            % Smoothing
            estimates = structfun(funMavg,estimates,'UniformOutput',false);
        end
    end 
end

end