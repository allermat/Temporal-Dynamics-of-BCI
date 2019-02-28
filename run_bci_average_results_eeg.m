function run_bci_average_results_eeg()
% Averages bci model results

expStage = 'final';
trMethod = 'sample-wise-sm-avg';

saveDf = cd(DEC_2_setupdir(expStage,'anal_eeg'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

fieldsOfInt = {'bic','logLike','bestParamFm'};

for i = 1:size(subjList,1)
    
    subID = subjList{i};
    fileMatchStr = ['bci_simul_MVPA_',subID,'_[0-9]+.mat'];
    saveDf = cd(fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID),trMethod));
    fileList = cellstr(ls);
    matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
    
    if sum(matchID) == 0
        warning('No file, skipping this subject! ');
        cd(saveDf);
        continue;
    elseif sum(matchID) > 1
        warning('Multiple files, skipping this subject! ');
        cd(saveDf);
        continue;
    else
        fileName = fileList(matchID);
        fileName = fileName{:};
    end
    
    if i == 1
        bciRes = load(fileName);
        out.timePoints = bciRes.timePoints;
        out.parameterNames = bciRes.parameterNames;
        out.neurDataSelection = bciRes.neurDataSelection;
        out.models = bciRes.models;
    else
        bciRes = load(fileName,fieldsOfInt{:});
    end
    
    for j = 1:numel(fieldsOfInt)
        temp = bciRes.(fieldsOfInt{j});
        tt = cell(1,size(temp,2));
        for k = 1:size(temp,2)
            tt{k} = cat(1,temp{:,k});
        end
        indivData{i,j} = tt;
    end
        
    cd(saveDf);
end


for i = 1:numel(fieldsOfInt)
    if strcmp(fieldsOfInt{i},'bestParamFm')
        % Average best parameters
        temp = cat(1,indivData{:,i});
        [tt,tts,tte] = deal(cell(1,size(temp,2)));
        for j = 1:size(temp,2)
            tt{j} = mean(cat(3,temp{:,j}),3);
            tts{j} = std(cat(3,temp{:,j}),0,3);
            tte{j} = tts{j}./sqrt(size(temp,1));
        end
        out.(fieldsOfInt{i}) = tt;
        out.([fieldsOfInt{i},'_rel']) = tt;
        out.([fieldsOfInt{i},'_std']) = tts;
        out.([fieldsOfInt{i},'_err']) = tte;
    elseif strcmp(fieldsOfInt{i},'bic')
        % Sum over BICs
        temp = cat(1,indivData{:,i});
        tt = cell(1,size(temp,2));
        for j = 1:size(temp,2)
            tt{j} = sum(cat(3,temp{:,j}),3);
        end
        out.([fieldsOfInt{i},'_fx']) = tt;
    end
end

% Sum BICs across participants


% Relative BIC
fieldName = 'bicRelBci';
temp = cell2mat(out.bic_fx);
temp = temp-repmat(out.bic_fx{1},1,size(temp,2));
out.(fieldName) = mat2cell(temp,size(temp,1),ones(size(temp,2),1));

fieldName = 'bicRelSegV';
temp = cell2mat(out.bic_fx);
temp = temp-repmat(out.bic_fx{4},1,size(temp,2));
out.(fieldName) = mat2cell(temp,size(temp,1),ones(size(temp,2),1));

% Saving subject specific bci simulations
fprintf('\n\nSaving data...\n\n');
savePath = fullfile(DEC_2_setupdir(expStage,'anal_eeg_group_mvpa'),trMethod,...
    'bci_simul_MVPA_group.mat');
save(savePath,'-struct','out','-v7.3');

end

