function stats = getStats(obj,varName,varargin)
% Method for accessing stats for various variables
%
% USAGE:
%   stats = getStats(obj,varName)
%   stats = getStats(obj,varName,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%       varName (string): name of the variable to which the
%           statistics belong
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       smooth (logical): whether to smooth the time course of AV estimates
% OUTPUT:
%   stats (struct): stats as fields of the array. Each field
%       contains a time series of stats for a certain condition, or a
%       time x time matrix of stats depending on input settings. 

% Copyright(C) 2017, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validGenTimes = {'tr','tr_x_tr'};
validVarNames = {'genPerf','avWeights','avModelCorr'};
addRequired(p,'obj');
addRequired(p,'varName',@(x)any(validatestring(x,validVarNames)));
addParameter(p,'genTime','tr',@(x)any(validatestring(x,validGenTimes)));
addParameter(p,'smooth',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
parse(p,obj,varName,varargin{:});
obj = p.Results.obj;
varName = p.Results.varName;
genTime = p.Results.genTime;
smooth = p.Results.smooth;

switch varName
  case 'genPerf'
    fieldName = 'stat_perfEstimates';
  case 'avWeights'
    if smooth
        fieldName = 'stat_AVmodelWeights_s';
    else
        fieldName = 'stat_AVmodelWeights';
    end
  case 'avModelCorr'
    if smooth
        fieldName = 'stat_AVmodelCorrelations_s';
    else
        fieldName = 'stat_AVmodelCorrelations';
    end
end
% Checking if AV model weight stats are part of the dataset
if any(~ismember(fieldName,obj.who))
    stats = [];
    warning('mvpares:getStats:reqestedDataNotPresent',...
        'The dataset does not contain stats for this variable. ');
    return;
end

% Loading AV model weihgt stats
if strcmp(genTime,'tr_x_tr')
    if smooth
        stats = [];
        warning('mvpares:getStats:invalidSettings',...
            'Smoothing is not applied to time x time generalized stats!');
        return;
    end
    stats = obj.data.(fieldName);
    statFieldNames = fieldnames(stats);
    if all(cellfun(@isempty,regexp(statFieldNames,'(h|p|st)_tr_x_tr_.*','once')))
        stats = [];
        warning('mvpares:getStats:reqestedDataNotPresent',...
            ['The dataset does not contain stats for this variable ',...
            'generalized time x time.']);
        return;
    else
        idx = cellfun(@isempty,regexp(statFieldNames,'(h|p|st)_tr_x_tr_.*','once'));
        toBeRemoved = statFieldNames(idx);
        stats = rmfield(stats,toBeRemoved);
        genTimeStr = 'tr_x_tr';
    end
elseif strcmp(genTime,'tr')
    stats = obj.data.(fieldName);

    statFieldNames = fieldnames(stats);
    if all(cellfun(@isempty,regexp(statFieldNames,'(h|p|st)_tr_(?!x_tr).*','once')))
        % Try if there are stats for across time generalized data
        % and use them if there are. 
        if any(~cellfun(@isempty,regexp(statFieldNames,...
                                       '(h|p|st)_tr_x_tr_.*','once')))
            idx = cellfun(@isempty,regexp(statFieldNames,...
                                          '(h|p|st)_tr_x_tr_.*','once'));
            toBeRemoved = statFieldNames(idx);
            stats = rmfield(stats,toBeRemoved);
            stats = structfun(@diag,stats,'UniformOutput',false);
            genTimeStr = 'tr_x_tr';
        else
            stats = [];
            warning('mvpares:getStats:reqestedDataNotPresent',...
                    ['The dataset does not contain stats for this variable ',...
                     'timeseries generalized.']);
            return;
        end
    else
        idx = cellfun(@isempty,regexp(statFieldNames,'(h|p|st)_tr_(?!x_tr).*','once'));
        toBeRemoved = statFieldNames(idx);
        stats = rmfield(stats,toBeRemoved);
        genTimeStr = 'tr';
    end
end
% Removing the generalization time flag from the field names
statFieldNames = fieldnames(stats);
for i = 1:numel(statFieldNames)
    stats.(strrep(statFieldNames{i},[genTimeStr,'_'],'')) = stats.(statFieldNames{i});
end
stats = rmfield(stats,statFieldNames);
stats = orderfields(stats);

end

