function [info,infoCorr,featCorr] = generateCorrectedAvData(info,feat,condDef,corrMode,dataType)
% Method for generating corrected AV data

if strcmp(dataType,'single_trial')
    % Selecting A and V only trials for averaging
    condsAxorV = condDef.condition(isnan(condDef.locationAuditory) | ...
                                   isnan(condDef.locationVisual));
    temp = info(ismember(info.condition,condsAxorV),{'condition', ...
                        'hand','task','locA','locV','relV'});
    [~,idx] = unique(temp(:,{'condition','hand'}),'rows');
    infoAvg = temp(idx,:);
    featAvg = NaN(size(feat,1),size(feat,2),size(infoAvg,1));

    for i = 1:size(infoAvg,1)
        actTrials = info.condition == infoAvg.condition(i) & ...
            info.hand == infoAvg.hand(i);
        featAvg(:,:,i) = mean(feat(:,:,actTrials),3);        
    end

    % Only AV trials will be included in the final dataset
    isAVtrial = ismember(info.condition,condDef.condition(~isnan(condDef.locationAuditory) & ...
                                                      ~isnan(condDef.locationVisual)));
    infoCorr = info(isAVtrial,{'session','iTrialSession','condition', ...
                        'hand','task','locA','locV','relV','resp'});
    featCorr = feat(:,:,isAVtrial);
    for i = 1:size(infoCorr,1)
        actAudAvgIdx = infoAvg.locA == infoCorr.locA(i) & ...
            infoAvg.hand == infoCorr.hand(i);
        actVisAvgIdx = infoAvg.locV == infoCorr.locV(i) & ...
            infoAvg.hand == infoCorr.hand(i) & ...
            infoAvg.relV == infoCorr.relV(i);
        % Correction mode 1: AV - (A + V)
        if corrMode == 1
            featCorr(:,:,i) = featCorr(:,:,i)-featAvg(:,:,actAudAvgIdx);
        end
        % Correction mode 2: AV - V
        featCorr(:,:,i) = featCorr(:,:,i)-featAvg(:,:,actVisAvgIdx);
    end
    
elseif strcmp(dataType,'avg_potential')

    % Only AV trials will be included in the final dataset
    % isAVtrial = ismember(info.condition,condDef.condition(~isnan(condDef.locationAuditory) & ...
                                                      % ~isnan(condDef.locationVisual)));
    infoCorr = info(:,{'session','condition','task','locA', ...
                        'locV','relV','resp'});
    featCorr = feat;
    for i = 1:size(infoCorr,1)
        % Correcting only AV examples, I keep the unisensory
        % examples for compatibility with the sample-wise-sm-avg
        % dataset
        if ~isnan(infoCorr.locA(i)) && ~isnan(infoCorr.locV(i))
            actAudIdx = info.locA == infoCorr.locA(i) & ...
                isnan(info.locV) & ...
                info.session == infoCorr.session(i);
            actVisIdx = info.locV == infoCorr.locV(i) & ...
                info.relV == infoCorr.relV(i) & ...
                isnan(info.locA) & ...
                info.session == infoCorr.session(i);
            % Correction mode 1: AV - (A + V)
            if corrMode == 1
                featCorr(:,:,i) = featCorr(:,:,i)-feat(:,:,actAudIdx);
            end
            % Correction mode 2: AV - V
            featCorr(:,:,i) = featCorr(:,:,i)-feat(:,:,actVisIdx);
        end
    end
    
else
    error('mvpa:generateCorrectedAvData:invalidInput',...
          'The value for dataType is invalid.');
end

end

