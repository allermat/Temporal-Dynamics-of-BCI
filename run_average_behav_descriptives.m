function run_average_behav_descriptives()
% Averages bci model results

expStage = 'final';

saveDf = cd(DEC_2_setupdir(expStage,'anal_behav'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

for i = 1:size(subjList,1)
    
    subID = subjList{i};
    fileMatchStr = ['descriptives_BEHAV_',subID,'.mat'];
    saveDf = cd(DEC_2_setupdir(expStage,'anal_behav_sub',subID));
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
    
    temp = load(fileName);
    indivData{i} = temp.respDistribution;
    
    cd(saveDf);
end
respDistribution = indivData{1};
indivData  = cat(1,indivData{:});

temp = varfun(@mean,indivData,...
    'InputVariables','distribution','GroupingVariables',...
    {'condition'});
respDistribution.GroupCount = temp.GroupCount;
respDistribution.distribution = temp.mean_distribution;

% Saving subject specific bci simulations
fprintf('\n\nSaving data...\n\n');
savePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_group'),...
    'descriptives_BEHAV_group.mat');
save(savePath,'respDistribution','-v7.3');

end

