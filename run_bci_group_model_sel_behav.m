function run_bci_group_model_sel_behav(varargin)
% Group level model selection for behavioural data

defaultModelSel = {'bci','fus','taskRel'};
validModelSel = {'bci','bciSel','bciMatch','fus','taskRel'};
p = inputParser;
addParameter(p,'modelSel',defaultModelSel,@(x) all(ismember(x,validModelSel)));
parse(p,varargin{:});

modelSel = p.Results.modelSel;

expStage = 'final';

saveDf = cd(DEC_2_setupdir(expStage,'anal_behav'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

fieldsOfInt = {'bestParamFm','bic','negLogLike','R2'};
nSubj = size(subjList,1);

for i = 1:nSubj
    
    subID = subjList{i};
    fileMatchStr = ['bci_simul_BEHAV_',subID,'_[0-9]+.mat'];
    saveDf = cd(DEC_2_setupdir(expStage,'anal_behav_sub',subID));
    listings = dir;
    fileList = {listings.name}';
    fileDates = [listings.datenum]';
    matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
        
    if sum(matchID) == 0
        cd(saveDf);
        error('No file for subject %s!',subID);
    elseif sum(matchID) > 1
        warning('Multiple files, using the most recent!');
        [~,idx] = max(fileDates(matchID));
        temp = fileList(matchID);
        fileName = temp(idx);
        fileName = fileName{:};
    else
        fileName = fileList(matchID);
        fileName = fileName{:};
    end
    
    if i == 1
        bciRes = load(fileName);
        mdlIdx = ismember(bciRes.models,modelSel);
        out.parameterNames = bciRes.parameterNamesAll(mdlIdx);
        out.models = bciRes.models(mdlIdx);
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

% Preparing negative log likelihoods for BMS
temp = cat(1,indivData{:,ismember(fieldsOfInt,'bic')});
temp = cell2mat(temp);
bicMat = temp;
bicMat = bicMat(:,ismember(out.models,modelSel),:);
[alpha,exp_r,xp,pxp,bor] = deal(NaN(size(bicMat,3),size(bicMat,2)));

for i = 1:size(bicMat,3)
    [alpha(i,:),exp_r(i,:),xp(i,:),pxp(i,:),bor(i,:)] = spm_BMS(bicMat(:,:,i));
end
out.rfx.alpha = alpha;
out.rfx.exp_r = exp_r;
out.rfx.xp = xp;
out.rfx.pxp = pxp;
out.rfx.bor = bor;

% Computing group level mean of the best parameters
temp = cat(1,indivData{:,ismember(fieldsOfInt,'bestParamFm')});
temp = arrayfun(@(x) cat(1,temp{:,x}),1:size(temp,2), ...
                'UniformOutput',false);
paramCell = temp;
paramCell = paramCell(:,ismember(out.models,modelSel),:);
out.rfx.paramMeans = cellfun(@mean,paramCell,'UniformOutput',false);
out.rfx.paramSTDs = cellfun(@std,paramCell,'UniformOutput',false);
out.rfx.paramSEMs = cellfun(@(x) x./sqrt(nSubj),out.rfx.paramSTDs, ...
                            'UniformOutput',false);

% Computing group level mean R2
temp = cat(1,indivData{:,ismember(fieldsOfInt,'R2')});
temp = cell2mat(temp);
temp = temp(:,ismember(out.models,modelSel),:);
out.rfx.R2Means = (tanh(mean(atanh(sqrt(temp))))).^2;
out.rfx.R2STDs = (tanh(std(atanh(sqrt(temp))))).^2;
out.rfx.R2SEMs = out.rfx.R2STDs./sqrt(nSubj);

% Sum over BICs for fixed effects analysis
temp = cat(1,indivData{:,ismember(fieldsOfInt,'bic')});
temp = temp(:,ismember(out.models,modelSel),:);
tt = cell(1,size(temp,2));
for j = 1:size(temp,2)
    tt{j} = sum(cat(3,temp{:,j}),3);
end
out.ffx.bic = tt;

% Relative BIC
fieldName = 'bicRelBci';
temp = cell2mat(out.ffx.bic);
temp = temp-repmat(out.ffx.bic{1},1,size(temp,2));
out.ffx.(fieldName) = mat2cell(temp,size(temp,1),ones(size(temp,2),1));

% fieldName = 'bicReltaskRel';
% temp = cell2mat(out.ffx.bic);
% temp = temp-repmat(out.ffx.bic{3},1,size(temp,2));
% out.ffx.(fieldName) = mat2cell(temp,size(temp,1),ones(size(temp,2),1));

% Saving subject specific bci simulations
fprintf('\n\nSaving data...\n\n');
nMdlStr = ['_',num2str(numel(modelSel)),'mdl'];
saveFileName = ['bci_model_sel_group',nMdlStr, ...
                '_',datestr(now,'yymmddHHMMSS'),'.mat'];
savePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_group'),...
                    saveFileName);
save(savePath,'-struct','out','-v7.3');

end

