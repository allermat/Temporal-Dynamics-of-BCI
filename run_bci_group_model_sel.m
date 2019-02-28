function run_bci_group_model_sel(trMethod,varargin)
% Group model selection

validTrMethods = {'sample-wise-sm','sample-wise-sm-avg'};
defaultModelSel = {'bci','fus','segA','segV','taskRel'};
p = inputParser;
addRequired(p,'trMethod',@(x) ismember(x,validTrMethods));
addParameter(p,'smooth',false,@islogical);
addParameter(p,'specIDstr','',@ischar);
addParameter(p,'modelSel',defaultModelSel,@(x) all(ismember(x,defaultModelSel)));
parse(p,trMethod,varargin{:});

trMethod = p.Results.trMethod;
smooth = p.Results.smooth;
specIDstr = p.Results.specIDstr;
modelSel = p.Results.modelSel;

expStage = 'final';

saveDf = cd(DEC_2_setupdir(expStage,'anal_eeg'));
subjList = cellstr(ls);
subjList = regexp(subjList,'[0-9]{3}','match','once');
subjList = subjList(~ismember(subjList,{''}));
cd(saveDf);

fieldsOfInt = {'bic','negLogLike','R2','plVar'};

nSubj = size(subjList,1);
if ~strcmp(specIDstr,'') && ~strcmp(specIDstr(1),'_')
    error('specIDstr must start with ''_''');
end

for i = 1:nSubj
    
    subID = subjList{i};
    fileMatchStr = ['bci_simul_MVPA_',subID,specIDstr,'_[0-9]+.mat'];
    saveDf = cd(fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID),trMethod,'control'));
    listings = dir;
    fileList = {listings.name}';
    fileDates = [listings.datenum]';
    matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
    
    if sum(matchID) == 0
        cd(saveDf);
        error('No file, for subject %s! ',subID);
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

% Checking temporal smoothing


% Preparing bic values for BMS
temp = cat(1,indivData{:,ismember(fieldsOfInt,'bic')});
if smooth
    % window of moving average, 4 samples equals 20 ms
    mavgWinSample = 4;
    % Kernel for smoothing trials with the given time wintow
    kern = ones(mavgWinSample,1)./mavgWinSample;
    % Moving average function
    funMavg = @(x) conv(x,kern,'same');
    temp = cellfun(funMavg,temp,'UniformOutput',false);
end
temp = cellfun(@(x) permute(x,[3,2,1]),temp,'UniformOutput',false);
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
out.modelSel = modelSel;

% Preparing R2s
temp = cat(1,indivData{:,ismember(fieldsOfInt,'R2')});
temp = cellfun(@(x) permute(x,[3,2,1]),temp,'UniformOutput',false);
temp = cell2mat(temp);
R2mat = temp;
R2mat = R2mat(:,ismember(out.models,modelSel),:);
out.rfx.R2Means =  squeeze((tanh(mean(atanh(sqrt(R2mat))))).^2)';
out.rfx.R2STDs = squeeze((tanh(std(atanh(sqrt(R2mat))))).^2)';
out.rfx.R2SEMs = out.rfx.R2STDs./sqrt(nSubj);

% Computing group level mean of the log-likelihoods
temp = cat(1,indivData{:,ismember(fieldsOfInt,'negLogLike')});
temp = cellfun(@(x) permute(x,[3,2,1]),temp,'UniformOutput',false);
temp = cell2mat(temp);
logLikemat = -temp;
logLikemat = logLikemat(:,ismember(out.models,modelSel),:);
out.rfx.logLikeMeans = squeeze(mean(logLikemat))';
out.rfx.logLikeSTDs = squeeze(std(logLikemat))';
out.rfx.logLikeSEMs = out.rfx.logLikeSTDs./sqrt(nSubj);

% Computing group level mean variance of the decoded labels
temp = cat(1,indivData{:,ismember(fieldsOfInt,'plVar')});
temp = cellfun(@transpose,temp,'UniformOutput',false);
plVarmat = cell2mat(temp);
out.rfx.plVarMean = mean(plVarmat)';
out.rfx.plVarSTD = std(plVarmat)';
out.rfx.plVarSEM = out.rfx.plVarSTD./sqrt(nSubj);

% Sum over BICs for fixed effects analysis
temp = cat(1,indivData{:,ismember(fieldsOfInt,'bic')});
temp = temp(:,ismember(out.models,modelSel));
tt = cell(1,size(temp,2));
for j = 1:size(temp,2)
    tt{j} = sum(cat(3,temp{:,j}),3);
end
out.ffx.bic = tt;

% Relative BIC
fieldName = 'bicRelBci';
temp = cell2mat(out.ffx.bic);
temp = temp-repmat(out.ffx.bic{ismember(modelSel,'bci')},1,size(temp,2));
out.ffx.(fieldName) = mat2cell(temp,size(temp,1),ones(size(temp,2),1));

fieldName = 'bicRelSegV';
temp = cell2mat(out.ffx.bic);
temp = temp-repmat(out.ffx.bic{ismember(modelSel,'segV')},1,size(temp,2));
out.ffx.(fieldName) = mat2cell(temp,size(temp,1),ones(size(temp,2),1));

% Saving subject specific bci simulations
fprintf('\n\nSaving data...\n\n');
if smooth, smStr = '_sm'; else smStr = ''; end
nMdlStr = ['_',num2str(numel(modelSel)),'mdl'];
saveFileName = ['bci_model_sel_MVPA_group',specIDstr,smStr,nMdlStr, ...
                '_',datestr(now,'yymmddHHMMSS'),'.mat'];
savePath = fullfile(DEC_2_setupdir(expStage,'anal_eeg_group_mvpa'),trMethod,'control',...
                    saveFileName);
save(savePath,'-struct','out','-v7.3');

end

