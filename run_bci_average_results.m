function run_bci_average_results()
% Averages bci model results

expStage = 'final';

saveDf = cd(DEC_2_setupdir(expStage,'anal_behav'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

for i = 1:size(subjList,1)
    
    subID = subjList{i};
    fileMatchStr = ['bci_simul_BEHAV_',subID,'_[0-9]+.mat'];
    saveDf = cd(DEC_2_setupdir(expStage,'anal_behav_sub',subID));
    fileList = cellstr(ls);
    matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
    
    if sum(matchID) == 0
        cd(saveDf);
        error('No file, skipping this subject! ');
    elseif sum(matchID) > 1
        cd(saveDf);
        error('Multiple files, skipping this subject! ');
    else
        fileName = fileList(matchID);
        fileName = fileName{:};
    end
    
    temp = load(fileName);
    
    indivData{i} = temp.mdlEval{1};
    
    cd(saveDf);
end
indivData = indivData';
fieldsOfInterest = {'freq_predV','freq_predA'};

bciSimulations = struct('conditions',{},'relV',{},'freq_predV',[],...
    'freq_predV_std',[],'freq_predV_sem',[],'freq_predA',[],...
    'freq_predA_std',[],'freq_predA_sem',[]);
bciSimulations(size(indivData{1},2)).conditions = {};

for i = 1:size(indivData{1},2)
    bciSimulations(i).conditions = indivData{1}(i).conditions;
    bciSimulations(i).relV = indivData{1}(i).relV;
    for j = 1:numel(fieldsOfInterest)
        extractFields = @(x) x(i).(fieldsOfInterest{j});
        temp = cell2mat(cellfun(extractFields,indivData,'UniformOutput',false));
        bciSimulations(i).(fieldsOfInterest{j}) = mean(temp);
        bciSimulations(i).([fieldsOfInterest{j},'_std']) = std(temp);
        bciSimulations(i).([fieldsOfInterest{j},'_sem']) = ...
            bciSimulations(i).([fieldsOfInterest{j},'_std'])./sqrt(size(temp,1));
    end
end


% Saving subject specific bci simulations
fprintf('\n\nSaving data...\n\n');
savePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_group'),...
    ['bci_simul_BEHAV_group','_',datestr(now,'yymmddHHMMSS'),'.mat']);
save(savePath,'bciSimulations','-v7.3');

end

