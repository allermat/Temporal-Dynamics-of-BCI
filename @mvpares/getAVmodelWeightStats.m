function stats = getAVmodelWeightStats(obj,varargin)
% Method for accessing AV model weights
%
% USAGE:
%   weigths = getAVmodelWeightStats(obj)
%   weigths = getAVmodelWeightStats(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       smooth (logical): whether to smooth the time course of AV estimates
% OUTPUT:
%   stats (struct): stats as fields of the array. Each field
%       contains a time series of weights for a certain condition, or a
%       time x time matrix of weights depending on input settings. 

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

% Checking if AV model weight stats are part of the dataset
if any(~ismember({'stat_AVmodelWeights','stat_AVmodelWeights_s'},obj.who))
    stats = [];
    warning('mvpares:getAVmodelWeightStats:reqestedDataNotPresent',...
        'The dataset does not contain AV model weight stats. ');
    return;
end

% Loading AV model weihgt stats
if strcmp(genTime,'tr_x_tr')
    if smooth
        stats = [];
        warning('mvpares:getAVmodelWeightStats:invalidSettings',...
            'Smoothing is not applied to time x time generalized weights!');
        return;
    end
    stats = obj.data.stat_AVmodelWeights;
    fieldNames = fieldnames(stats);
    if all(cellfun(@isempty,regexp(fieldNames,'(h|p|st)_tr_x_tr_.*','once')))
        stats = [];
        warning('mvpares:getAVmodelWeightStats:reqestedDataNotPresent',...
            ['The dataset does not contain AV model weight stats for ',...
            'time x time generalized weights.']);
        return;
    else
        idx = cellfun(@isempty,regexp(fieldNames,'(h|p|st)_tr_x_tr_.*','once'));
        toBeRemoved = fieldNames(idx);
        stats = rmfield(stats,toBeRemoved);
    end
elseif strcmp(genTime,'tr')
    if smooth
        stats = obj.data.stat_AVmodelWeights_s;
    else
        stats = obj.data.stat_AVmodelWeights;
    end
    fieldNames = fieldnames(stats);
    if all(cellfun(@isempty,regexp(fieldNames,'(h|p|st)_tr_(?!x_tr).*','once')))
        stats = [];
        warning('mvpares:getAVmodelWeightStats:reqestedDataNotPresent',...
            ['The dataset does not contain AV model weight stats for ',...
            'timeseries generalized weights.']);
        return;
    else
        idx = cellfun(@isempty,regexp(fieldNames,'(h|p|st)_tr_(?!x_tr).*','once'));
        toBeRemoved = fieldNames(idx);
        stats = rmfield(stats,toBeRemoved);
    end
end
% Removing the generalization time flag from the field names
fieldNames = fieldnames(stats);
for i = 1:numel(fieldNames)
    stats.(strrep(fieldNames{i},[genTime,'_'],'')) = stats.(fieldNames{i});
end
stats = rmfield(stats,fieldNames);
stats = orderfields(stats);

end

