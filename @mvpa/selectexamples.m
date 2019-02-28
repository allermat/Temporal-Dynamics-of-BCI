function isExample = selectexamples(inputConds,condDef,dataInfo)
% Selects examples matching the given input condition
% 
% INPUT:
%   cond: hyphen separated string of the input conditions
%   condDef: table of stimulus condition definitions
%   dataInfo: table of information about the examples
% 
% OUTPUT:
%   isExample: N x 1 logical vector indicating which examples correspond to
%       the given input conditions. N = number of examples in dataInfo.
%

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

%% Parsing input
p = inputParser;

validConds = {'A','V','AV','AxorV','pseudoAV','c','ci','i','a','v','av','hrel','lrel',...
    'lh','rh','la','ha'};

addRequired(p,'inputConds',@(x)all(ismember(regexp(x,'-','split'),validConds)));
addRequired(p,'condDef',@(x)validateattributes(x,{'table'},{'nonempty'}));
addRequired(p,'dataInfo',@(x)validateattributes(x,{'table'},{'nonempty'}));

% Parsing inputs
parse(p,inputConds,condDef,dataInfo);

% Assigning inputs to variables
inputConds = p.Results.inputConds;
condDef = p.Results.condDef;
dataInfo = p.Results.dataInfo;

%% Main
% Splitting the input conditions to parts
condsCell = regexp(inputConds,'-','split');
% Making sure, that only one level of a condition category is specified at 
% a time
moreThanOnePerCat = sum(ismember({'A','V','AV','AxorV','pseudoAV'},condsCell)) > 1 || ...
                    sum(ismember({'c','i','ci'},condsCell)) > 1 || ...
                    sum(ismember({'a','v','av'},condsCell)) > 1 || ...
                    sum(ismember({'lrel','hrel'},condsCell)) > 1 || ...
                    sum(ismember({'lh','rh'},condsCell)) > 1 || ...
                    sum(ismember({'la','ha'},condsCell)) > 1;
if moreThanOnePerCat
    error('mvpa:selectexamples:ambiguousCondition',...
        'Only one level of a condition category can be specified at a time');
end
% Making sure, that the modality condition is at the beginning
if any(ismember({'A','V','AV','AxorV','pseudoAV'},condsCell))
    temp = condsCell(ismember(condsCell,{'A','V','AV','AxorV','pseudoAV'}));
    condsCell(ismember(condsCell,{'A','V','AV','AxorV','pseudoAV'})) = [];
    condsCell = [temp,condsCell];
end
% Making sure, that the conditions which are assessed in situ (e.g. low vs.
% high prestimulus alpha) are at the end of the condition list, so they are
% evaluated last. 
if any(ismember({'la','ha'},condsCell))
    temp = condsCell(ismember(condsCell,{'la','ha'}));
    condsCell(ismember(condsCell,{'la','ha'})) = [];
    condsCell = [condsCell,temp];
end

% Flag indicating if the modality condition is pseudoAV. This is
% necessary as the congruency and response modality conditions are
% not evaluated here if it is pseudoAV. 
isPseudoAV = false;
isAxorV = false;

% Logical array keeping track of the selected examples based on the conditions
isExample = true(size(dataInfo,1),1);

for i = 1:size(condsCell,2)
    
    switch condsCell{i}
        
      case 'A'
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(isnan(condDef.locationVisual)));
      case 'V'
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(isnan(condDef.locationAuditory)));
      case 'AV'
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(~isnan(condDef.locationAuditory) & ...
                                                          ~isnan(condDef.locationVisual)));
      case 'AxorV'
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(isnan(condDef.locationAuditory) | ...
                                                          isnan(condDef.locationVisual)));
        isAxorV = true;
      case 'pseudoAV'
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(isnan(condDef.locationAuditory) | ...
                                                          isnan(condDef.locationVisual)));
        isPseudoAV = true;
      case 'c'
        if isPseudoAV, continue; end
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(condDef.locationAuditory == condDef.locationVisual));
      case 'i'
        if isPseudoAV, continue; end
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(condDef.locationAuditory ~= condDef.locationVisual));
      case 'ci'
        if isPseudoAV, continue; end
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(condDef.locationAuditory == condDef.locationVisual | ...
                                                          condDef.locationAuditory ~= condDef.locationVisual));
      case 'a'
        if isPseudoAV, continue; end
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(condDef.task == 1));
      case 'v'
        if isPseudoAV, continue; end
        isExample = isExample & ismember(dataInfo.condition,...
                                         condDef.condition(condDef.task == 2));
      case 'av'
        if isPseudoAV, continue; end
        isExample = isExample & ismember(dataInfo.condition,...
                                         [condDef.condition(condDef.task == 1);condDef.condition(condDef.task == 2)]);
      case 'lrel'
        if isAxorV || isPseudoAV
            % In these cases we also include A trials as well even
            % if they don't ave visual reliability information
            isExample = isExample & ...
                (isnan(dataInfo.locV) | dataInfo.relV == 12);
        else
            isExample = isExample & dataInfo.relV == 12;
        end
      case 'hrel'
        if isAxorV || isPseudoAV
            % In these cases we also include A trials as well even
            % if they don't ave visual reliability information
            isExample = isExample & ...
                (isnan(dataInfo.locV) | dataInfo.relV == 2);
        else
            isExample = isExample & dataInfo.relV == 2;
        end
      case 'lh'
        isExample = isExample & dataInfo.hand == 1;
      case 'rh'
        isExample = isExample & dataInfo.hand == 2;
      case 'la'
        % Median-splitting the trials based on the specified conditions
        % so far
        m = nanmedian(dataInfo.psAlphaPowOcc(isExample));
        isExample = isExample & dataInfo.psAlphaPowOcc <= m;
      case 'ha'
        m = nanmedian(dataInfo.psAlphaPowOcc(isExample));
        isExample = isExample & dataInfo.psAlphaPowOcc > m;
      otherwise
        error('Unrecognized condition');     
    end
    
end

end