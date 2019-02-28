function weights = getAVmodelWeights(obj,varargin)
% Method for accessing AV model weights
%
% USAGE:
%   weigths = getAVmodelWeights(obj)
%   weigths = getAVmodelWeights(obj,'Name',Value)
% INPUT:
%   Required:
%       obj (object): mvpares object
%   'Name'-Value arguments:
%       genTime (string): indicating the generalization time. Possible 
%           values: 'tr','tr_x_tr' for trainin time and training time by 
%           training time respectively. Default: 'tr'. 
%       smooth (logical): whether to smooth the time course of AV estimates
%       unit (string): the unit of the output, 'deg' - degreees, 
%           'rad' - radians, default: 'rad'
% OUTPUT:
%   weights (struct): weights as fields of the array. Each field
%       contains a time series of weights for a certain condition, or a
%       time x time matrix of weights depending on input settings. 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validGenTimes = {'tr','tr_x_tr'};
validUnits = {'deg','rad'};
addRequired(p,'obj');
addParameter(p,'genTime','tr',@(x)any(validatestring(x,validGenTimes)));
addParameter(p,'smooth',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'unit','rad',@(x)any(validatestring(x,validUnits)));
parse(p,obj,varargin{:});
obj = p.Results.obj;
genTime = p.Results.genTime;
smooth = p.Results.smooth;
unit = p.Results.unit;

if strcmp(obj.level,'group')
    % Group level AV weigths are stored in the dataset
    
    % Checking if AV model weigths are part of the dataset
    if any(~ismember({'gen_AVmodelWeights','gen_AVmodelWeights_s'},obj.who))
        weights = [];
        warning('mvpares:getAVmodelWeights:reqestedDataNotPresent',...
            'The dataset does not contain AV model weights. ');
        return;
    end
    
    if strcmp(genTime,'tr_x_tr')
        % Check if genTime is valid given the result dataset
        if isvector(obj.getGenTimePoints)
            weights = [];
            warning('mvpares:getAVmodelWeights:reqestedDataNotPresent',...
                'The dataset is not generalized time x time! ');
            return;
        end
        if smooth
            weights = [];
            warning('mvpares:getAVmodelWeights:invalidSettings',...
                ['Can''t smooth if time x time generalized estimates ',...
                'are required! ']);
            return;
        end
        weights = obj.data.gen_AVmodelWeights;
    elseif strcmp(genTime,'tr')
        if smooth
            weights = obj.data.gen_AVmodelWeights_s;
        else
            weights = obj.data.gen_AVmodelWeights;
            % Extracting the diagonal if just the training timepoints are 
            % needed and the data is time x time generalized
            if ~isvector(obj.getGenTimePoints)
                weights = structfun(@diag,weights,'UniformOutput',false);
            end
        end
    end
        
else
    
    % Subject level AV weigths are computed on-line from the AV estimates
    estimates = obj.getAVmodelEstimates('genTime',genTime,'smooth',smooth);
    
    if isempty(estimates)
        warning('mvpares:getAVmodelWeights:reqestedDataNotPresent',...
            'Couldn''t find AV model estimates in the dataset, returning.');
        weights = [];
        return;
    end
    
    estimateNames = fieldnames(estimates);
    % Get the name of each present condition
    temp = estimateNames(~cellfun(@isempty,regexp(estimateNames,'A_.*_beta')));
    condNamesAll = strrep(strrep(temp,'A_',''),'_beta','');
    % All Wavs are computed from the 3-way interaction betas, so I
    % just select those betas
    idx = ~cellfun(@isempty,regexp(estimateNames,'._._._._.*'));
    estimateNames = estimateNames(idx);
    estimatesCell = struct2cell(estimates);
    estimatesCell = estimatesCell(idx);
    condNames3way = condNamesAll(~cellfun(@isempty,regexp(condNamesAll,'._._.')));
    
    % Separating A and V betas, matched by alphabetical order
    aBetas = estimatesCell(~cellfun(@isempty,regexp(estimateNames,'A_.*_beta')));
    vBetas = estimatesCell(~cellfun(@isempty,regexp(estimateNames,'V_.*_beta')));
    aCI = estimatesCell(~cellfun(@isempty,regexp(estimateNames,'A_.*_ci')));
    vCI = estimatesCell(~cellfun(@isempty,regexp(estimateNames,'V_.*_ci')));
        
    % Computing Wavs from the 3-way interaction term betas
    funb2w = @(v,a) atan2(v,a);
    wavs3way = cellfun(funb2w,vBetas,aBetas,'UniformOutput',false);
    aCIhigh = cellfun(@plus,aBetas,aCI,'UniformOutput',false);
    vCIhigh = cellfun(@plus,vBetas,vCI,'UniformOutput',false);
    wavsCIhigh = cellfun(funb2w,vCIhigh,aCIhigh,'UniformOutput',false);
    wavsCI3way = cellfun(@minus,wavsCIhigh,wavs3way,'UniformOutput',false);
    
    % Pre-allocating arrays for output data
    wavNamesAll = strcat(condNamesAll,repmat({'_wav'},size(condNamesAll)));
    ciNamesAll = strcat(condNamesAll,repmat({'_ci'},size(condNamesAll)));
    [wavsAll,wavsCIAll] = deal(cell(size(wavNamesAll)));
    
    % Computing Wavs for 2-way interaction and main effect terms    
    % Find 2-way interaction and main effect terms
    not3wayCond = ~ismember(condNamesAll,condNames3way);
    dimToCat = sum(size(estimatesCell{1}) > 1)+1;
    for i = 1:numel(not3wayCond)
        if not3wayCond(i)
            % Selecting 3 way interaction terms corresponding to
            % the 2-way or main effect term at hand
            actCond = strrep(condNamesAll{i},'_','');
            matchStr = ['[',actCond,']'];
            idx = cellfun(@(x) numel(x) == numel(actCond), ....
                          regexp(condNames3way,matchStr));
            
            % Averaging over the selected terms
            temp = cat(dimToCat,wavs3way{idx});
            wavsAll{i} = circ_mean(temp,[],dimToCat);
            temp = cat(dimToCat,wavsCI3way{idx});
            wavsCIAll{i} = circ_mean(temp,[],dimToCat);
        else
            wavsAll{i} = wavs3way{ismember(condNames3way, ...
                                           condNamesAll(i))};
            wavsCIAll{i} = wavsCI3way{ismember(condNames3way, ...
                                           condNamesAll(i))};
        end
    end
    
    weights = cell2struct(cat(1,wavsAll,wavsCIAll),cat(2,wavNamesAll,ciNamesAll));
    
end

% Converting to degrees if necessary
if strcmp(unit,'deg')
    weights = structfun(@degrees,weights,'UniformOutput',false);
end

end

