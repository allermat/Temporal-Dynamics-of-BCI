function mvparesObj = addGenPerfEstimates(mvparesObj)
% Method for computing generalization performance estimates
% 
% USAGE:
%   mvparesObj = addGenPerfEstimates(mvparesObj)
% INPUT:
%   mvparesObj (object): mvpares object
% OUTPUT:
%   mvparesObj (object): mvpares object with genaralization performance 
%       estimates

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
addRequired(p,'mvparesObj',@(x) isa(x,'mvpares') && x.isvalid);
parse(p,mvparesObj);
mvparesObj = p.Results.mvparesObj;

nFolds = mvparesObj.getNcvFolds;
sGenTime = mvparesObj.getSizeGenTime;
info = mvparesObj.getInfo;

if sGenTime(2) > 1
    genTimeStr = 'tr_x_tr';
else
    genTimeStr = 'tr';
end
% Opening parallel pool.
% if there is no parallel pool running, open one.
currPool = gcp('nocreate');
if isempty(currPool)
    parpool('local');
end

% Computing prediction performance pooled over folds
predLabels = mvparesObj.getPredLabels('genTime',genTimeStr);
trueLabels = mvparesObj.getTrueLabels;
if isempty(predLabels) || isempty(trueLabels)
    warning('mvpa:addGenPerfEstimates:requestedDataNotPresent',...
        ['Couldn''t extract true or predicted labels, performance estimates ',...
        'can''t be estimated, returning.']);
    return;
end

% Starting the timer and printing details
cStart = clock;
ds = datestr(now);
printFnTitle(80,'addGenPerfEstimates',ds)
fprintf('Estimating... \n');

if ismember(info.svm_type,{'cc','nc'})
    acc = mvpa.compPerfEstimates(trueLabels,predLabels,'classification');
    temp = struct();
    temp.acc_poolFolds = acc;
    
    % Computing prediction performance for each fold separately
    acc = NaN(sGenTime(1),sGenTime(2),nFolds);
    for i = 1:nFolds
        predLabels = mvparesObj.getPredLabels('cvFold',i,'genTime',genTimeStr);
        trueLabels = mvparesObj.getTrueLabels('cvFold',i);
        [acc(:,:,i)] = mvpa.compPerfEstimates(trueLabels,predLabels,'classification');
    end
    temp.acc_perFolds = acc;
    % Averaging over folds
    temp.acc_avgFolds = mean(acc,3);
else
    [int,b,r2,r] = mvpa.compPerfEstimates(trueLabels,predLabels,'regression');
    temp = struct();
    temp.int_poolFolds = int;
    temp.b_poolFolds = b;
    temp.r2_poolFolds = r2;
    temp.r_poolFolds = r;
    temp.rf_poolFolds = atanh(r);
    
    % Computing prediction performance for each fold separately
    [int,b,r2,r] = deal(NaN(sGenTime(1),sGenTime(2),nFolds));
    for i = 1:nFolds
        predLabels = mvparesObj.getPredLabels('cvFold',i,'genTime',genTimeStr);
        trueLabels = mvparesObj.getTrueLabels('cvFold',i);
        [int(:,:,i),b(:,:,i),r2(:,:,i),r(:,:,i)] = mvpa.compPerfEstimates(trueLabels,predLabels,'regression');
    end
    temp.int_perFolds = int;
    temp.b_perFolds = b;
    temp.r2_perFolds = r2;
    temp.r_perFolds = r;
    temp.rf_perFolds = atanh(r);
    
    % Averaging over folds
    int = mean(int,3);
    b = mean(b,3);
    % Fisher transform first, then average when averaging correlation
    % coefficients. 
    r2 = squeeze((tanh(mean(atanh(sqrt(r2)),3))).^2);
    r = squeeze(tanh(mean(atanh(r),3)));
    temp.int_avgFolds = int;
    temp.b_avgFolds = b;
    temp.r2_avgFolds = r2;
    temp.r_avgFolds = r;
    temp.rf_avgFolds = squeeze(mean(temp.rf_perFolds,3));
end

% Setting the dataset object to writable if it is not
if ~mvparesObj.writable, mvparesObj.setWritable(true); end

fieldName = 'gen_perfEstimates';
if ismember(fieldName,fieldnames(mvparesObj.data))
    warning('mvpa:addGenPerfEstimates:overwriteField',...
        ['The field ''%s'' already exists in the ',...
        'mvpa result dataset, it will be overwritten. '],fieldName);
end
mvparesObj.data.(fieldName) = temp;
mvparesObj.setWritable(false);

% Finishing timer and printing elapsed time
fprintf('Estimation elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

end