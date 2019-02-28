function [excProb,excProbNorm,indivBCIdata,indivBCIdataNorm] = comp_exceedance_prob_bci(varargin)

% Parsing input, checking matlab
p = inputParser;

validPerfIndices = {'r','r2','b','acc'};
validCvEval = {'avgFolds','poolFolds'};

addParameter(p,'perfIdx','r',@(x) any(validatestring(x,validPerfIndices)));
addParameter(p,'cvEval','avgFolds',@(x) any(validatestring(x,validCvEval)));

parse(p,varargin{:});

perfIdx = p.Results.perfIdx;
cvEval = p.Results.cvEval;

% Loading data
expStage = 'final';
trainMethod = 'sample-wise-sm-avg';

addpath(DEC_2_setupdir(expStage,'anal_scripts'));

saveDf = cd(DEC_2_setupdir(expStage,'anal_eeg'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

estIdxStr = [perfIdx,'_',cvEval];
matchStr = 'er_tr-AV-ci-[av]+_([a-zA-Z]+)_-100-5-1000_gen-AV-ci-[av]+_.*';
estimateStr = {'sAhatInd','sVhatInd','sHatComm','sAhatANDsVhat'};
[indivBCIdata,indivBCIdataNorm] = deal(cell(1,size(subjList,1)));

for iSub = 1:size(subjList,1)
    
    subID = subjList{iSub};
    saveDf = cd(fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID),trainMethod));
    
    % Choosing appropriate files
    fileList = cellstr(ls('*.mat'));
    matchIdx = ~cellfun(@isempty,regexp(fileList,matchStr));
    fileList = fileList(matchIdx);
    % Finding tokens indicating different BCI estimates
    estimateTokens = regexp(fileList,matchStr,'tokens');
    estimateTokens = [estimateTokens{:}]';
    estimateTokens = [estimateTokens{:}]';
    % Re-ordering the list of files according to the desired frequency order
    [~,Idx] = ismember(estimateStr,estimateTokens);
    fileList = fileList(Idx);
        
    if numel(fileList) ~= numel(estimateStr)
        error('Missing file! ');
    end
    
    temp = [];
    
    for iFile = 1:numel(fileList)
        mvparesObj = mvpares(fileList{iFile});
        perfEstimates = getGenPerfEstimates(mvparesObj,'genTime','tr');
        temp(:,iFile) = perfEstimates.(estIdxStr);
    end
    indivBCIdata{iSub} = temp;
    m = max(temp);
    % Adding a very small number just to avoid having exactly 1 at the
    % maximum
    m = m + 1.5*eps(0.5);
    indivBCIdataNorm{iSub} = temp./repmat(m,size(temp,1),1);
    cd(saveDf);
end

% Reshaping data to best suit the bootstrapping
temp = cat(3,indivBCIdata{:});
temp = mat2cell(temp,ones(size(temp,1),1),ones(size(temp,2),1),size(temp,3));
indivBCIdata = cellfun(@squeeze,temp,'UniformOutput',false);
temp = cat(3,indivBCIdataNorm{:});
temp = mat2cell(temp,ones(size(temp,1),1),ones(size(temp,2),1),size(temp,3));
indivBCIdataNorm = cellfun(@squeeze,temp,'UniformOutput',false);
% Number of bootstrap samples
nBoot = 10000;
% Function for determining the maximum location of an ordered series of 
% inputs
if any(ismember(perfIdx,{'r','r2'}))
    bootFun = @(x1,x2,x3,x4) tanh(mean(atanh([x1,x2,x3,x4]))) == max(tanh(mean(atanh([x1,x2,x3,x4]))));
elseif strcmp(perfIdx,'b')
    bootFun = @(x1,x2,x3,x4) mean([x1,x2,x3,x4]) == max(mean([x1,x2,x3,x4]));
else
    error('This performance index is not yet supported');
end
% Arrays for collecting exceedance probabilities
[excProb,excProbNorm] = deal(NaN(size(indivBCIdata)));

% Computing exceedance probability for raw accuracies
fprintf('Bootstrapping raw accuracies... \n');
parfor_progress(size(indivBCIdata,1));
for iSample = 1:size(indivBCIdata,1)
    actSample = indivBCIdata(iSample,:);
    bootStat = bootstrp(nBoot,bootFun,actSample{:},'Options',statset('UseParallel',true));
    excProb(iSample,:) = sum(bootStat)./nBoot;
    parfor_progress;
end
parfor_progress(0);

% Computing exceedance probability for normalized accuracies
fprintf('Bootstrapping normalized accuracies... \n');
parfor_progress(size(indivBCIdata,1));
for iSample = 1:size(indivBCIdataNorm,1)
    actSample = indivBCIdataNorm(iSample,:);
    bootStat = bootstrp(nBoot,bootFun,actSample{:},'Options',statset('UseParallel',true));
    excProbNorm(iSample,:) = sum(bootStat)./nBoot;
    parfor_progress;
end
parfor_progress(0);

colors = {[0,176,68]./255,[170,0,0]./255,[255,184,69]./255,[0,69,177]./255};
time = getTrTimePoints(mvparesObj)*1000;

figure();
subplot(1,2,1);
h = area(time,excProb);
for i = 1:numel(h)
    set(h(i),'FaceColor',colors{i});
end
ylim([0,1]);
title('Exceedance probability on raw betas');
xlabel('time (ms)');

subplot(1,2,2);
h = area(time,excProbNorm);
for i = 1:numel(h)
    set(h(i),'FaceColor',colors{i});
end
ylim([0,1]);
title('Exceedance probability on normalized betas');
xlabel('time (ms)');
legend(estimateStr);


rmpath(DEC_2_setupdir(expStage,'anal_eeg_scripts'));