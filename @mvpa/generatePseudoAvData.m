function [infoPseu,featPseu,groupingsPseu] = generatePseudoAvData(info,feat,groupings,trCond,condDef)

% Making sure we only use A and V only trials
if any(~isnan(info.locV) & ~isnan(info.locA))
    error('mvpa:generatePseudoAvData:invalidInput', ... 
          'Input data must only contain A only and V only conditions');
end

groups = unique(groupings(:,1,1));
info.group = groupings(:,1,1);

locAlevels = unique(info.locA(~isnan(info.locA)));
locVlevels = unique(info.locV(~isnan(info.locV)));

% Pre-allocating arrays
[featPseuTemp,infoPseuTemp,groupingsPseuTemp] = ...
    deal(cell(length(locAlevels),length(locVlevels),length(groups)));
% Making sure all output conditions have roughly equal number of
% examples
nPartsOut = max([numel(locAlevels),numel(locVlevels)]);
    function out = fun(x)
        % This makes sure that not always the numbers close to
        % nPartsOut drop out
        out = mod((1:size(x,1))'+randi(nPartsOut),nPartsOut);
        out(out == 0) = nPartsOut;
    end
% By calling the function separately for each condition I make
% sure, that each part is going to contain roughly equal number of
% high and low reliable visual trials in case both are present
temp = varfun(@fun,info,'InputVariables','condition','GroupingVariables',...
    {'condition','group'},'OutputFormat','cell');
info.part = cell2mat(temp);

for i = 1:numel(groups)
    vPartIdx = 1;
    for j = 1:numel(locAlevels)
        aPartIdx = 1;
        for k = 1:numel(locVlevels)
            trialIdxA = find(info.group == groups(i) & ...
                             info.locA == locAlevels(j) & ...
                             info.part == aPartIdx);
            trialIdxV = find(info.group == groups(i) & ...
                             info.locV == locVlevels(k) & ...
                             info.part == vPartIdx);
            
            % Resampling trials if necessary
            if size(trialIdxA,1) < size(trialIdxV,1)
                trialIdxV = trialIdxV(randsample(size(trialIdxV,1), ...
                                                 size(trialIdxA,1)));
            else % if size(trialIdxV,1) < size(trialIdxA,1)
                trialIdxA = trialIdxA(randsample(size(trialIdxA,1), ...
                                                 size(trialIdxV,1)));
            end
            
            featPseuTemp{j,k,i} = sum(cat(4,feat(:,:,trialIdxA),feat(:,:,trialIdxV)),4);
            tempInfo = info(trialIdxV,{'condition','hand','task','locA','locV','relV','resp'});
            tempInfo.locA = info.locA(trialIdxA);
            infoPseuTemp{j,k,i} = tempInfo;
            groupingsPseuTemp{j,k,i} = repmat(groups(i),size(tempInfo,1),1);
            
            if any(isnan(tempInfo.locA))
                error('NaN');
            end
            
            aPartIdx = aPartIdx + 1;
        end
        
        vPartIdx = vPartIdx + 1;
    end
end

infoPseu = cat(1,infoPseuTemp{:});
featPseu = cat(3,featPseuTemp{:});
groupingsPseu = cat(1,groupingsPseuTemp{:});

infoPseu.task = NaN(size(infoPseu.task));
infoPseu.hand = NaN(size(infoPseu.hand));
infoPseu.resp = NaN(size(infoPseu.resp));
for i = 1:size(infoPseu,1)
    % Assigning condition labels to the pseudoAV trials. I use the
    % condition labels corresponding to conditions with the given
    % locA, locV and relV and keep  the task as auditory
    % localization as in the generalization this info 
    % will not be specified. This is a workaround, might 
    % require a proper solution. 
    infoPseu.condition(i) = condDef.condition(...
        condDef.task == 1 & ...
        condDef.locationAuditory == infoPseu.locA(i) & ...
        condDef.locationVisual == infoPseu.locV(i) & ...
        condDef.reliabilityVisual == infoPseu.relV(i));
end
infoPseu.example = (1:size(infoPseu,1))';
infoPseu = infoPseu(:,[end,1:end-1]);

% Selecting trials based on training condition (congruency)
trCondsCell = regexp(trCond,'-','split');
if ismember('c',trCondsCell)
    isExample = ismember(infoPseu.condition, ...
                         condDef.condition(condDef.locationAuditory == ...
                                           condDef.locationVisual));
elseif ismember('i',trCondsCell)
        isExample = ismember(infoPseu.condition, ...
                         condDef.condition(condDef.locationAuditory ~= ...
                                           condDef.locationVisual));
else
    isExample = true(size(infoPseu,1),1);
end

infoPseu = infoPseu(isExample,:);
featPseu = featPseu(:,:,isExample);
groupingsPseu = groupingsPseu(isExample);
groupingsPseu = repmat(groupingsPseu,1,size(groupings,2), ...
                       size(groupings,3));

end

