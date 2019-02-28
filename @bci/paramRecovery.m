function paramRecovery(I)
% Method for running parameter recovery for a specified model and parameters

% Parsing input
p = inputParser;

validModels = {'bci','fus','segA','segV','taskRel'};
validResp = {'cont','discr'};
validDecisionFun = {'mdlAvg','mdlSel','probMatch'};

addParameter(p,'model','bci',@(x) any(validatestring(x,validModels)));
addParameter(p,'trueParam',struct,@(x) validateattributes(x,{'struct'},{'nonempty'}));
addParameter(p,'nCond',40,@(x) validateattributes(x,{'numeric'},...
    {'integer','positive'}));
addParameter(p,'nTrialPerCond',100,@(x) validateattributes(x,{'numeric'},...
    {'integer','positive'}));
addParameter(p,'resp','cont',@(x) any(validatestring(x,validResp)));
addParameter(p,'decisionFun','',@(x) any(validatestring(x,validDecisionFun)));
addParameter(p,'isServer',false,@(x) validateattributes(x,{'logical'},...
    {'scalar'}));

parse(p,I);

model = p.Results.model;
trueParam = p.Results.trueParam;
nCond = p.Results.nCond;
nTrialPerCond = p.Results.nTrialPerCond;
resp = p.Results.resp;
decisionFun = p.Results.decisionFun;
isServer = p.Results.isServer;

% Printing header for console output
ds = datestr(now);
printFnTitle(80,'run_bci_param_recovery',ds)
cStart = clock;

if strcmp(model,'fus')
    % Paramter settings for fake data generation
    cfg.model = model;
    cfg.param.sigP = trueParam.sigP;
    cfg.param.sigA = trueParam.sigA;
    cfg.param.sigV = trueParam.sigV;
    cfg.param.kW = trueParam.kW;
    cfg.param.muP = trueParam.muP;
    cfg.nCond = nCond;
    cfg.nTrialPerCond = nTrialPerCond;
    
    psDataTr = bci.generateFakeResponses(cfg);
    psDataTe = bci.generateFakeResponses(cfg);
    
    % Parameter settings for optimization
    n = 5;
    sigV = linspace(0.1,30,n);
    sigA = linspace(0.1,30,n);
    sigP = linspace(0.1,30,n);
    kW = 1;
    
    % Generating parameter combinations for full model
    parameterNames = {'sigP','sigA','sigV','kW'};
    gridVectors = {sigP,sigA,sigV,kW};
    fitFun = @bci.fitModelFus;
    trueParamVect = [trueParam.sigP,trueParam.sigA,trueParam.sigV,trueParam.kW];
    
elseif strcmp(model,'bci')
    % Paramter settings for fake data generation
    cfg.model = model;
    cfg.param.p_common = trueParam.p_common;
    cfg.param.sigP = trueParam.sigP;
    cfg.param.sigA = trueParam.sigA;
    cfg.param.sigV = trueParam.sigV;
    cfg.param.kW = trueParam.kW;
    cfg.param.muP = trueParam.muP;
    cfg.nCond = nCond;
    cfg.nTrialPerCond = nTrialPerCond;
    
    psDataTr = bci.generateFakeResponses(cfg);
    psDataTe = bci.generateFakeResponses(cfg);
    
    % Parameter settings for optimization
    n = 5;
    p_common = linspace(0,1,n);
    sigV = linspace(0.1,30,n);
    sigA = linspace(0.1,30,n);
    sigP = linspace(0.1,30,n);
    kW = 1;
    
    % Generating parameter combinations for full model
    parameterNames = {'p_common','sigP','sigA','sigV','kW'};
    gridVectors = {p_common,sigP,sigA,sigV,kW};
    fitFun = @(x1,x2,x3,x4) bci.fitModel(x1,x2,x3,x4,1);
    trueParamVect = [trueParam.p_common,trueParam.sigP,trueParam.sigA,...
        trueParam.sigV,trueParam.kW];
else
    error('bci:paramRecovery:missingModel','This model has not yet been implemented');
end

nParameters = numel(gridVectors);
% Full factorial expansion of the specified parameter vectors
coords = cell(1,nParameters);
[coords{:}] = ndgrid(gridVectors{:});
coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
% Matrix of all possible parameter combinations: nRows = number of
% combinations, columns correspond to parameters
paramCombinations = cat(2,coords{:});

% Pre-allocating array for data collection
logLike_all = NaN(size(paramCombinations,1),1);

fprintf('Performing gridsearch...\n');
if ~isServer, parfor_progress(size(paramCombinations,1)); end
parfor i = 1:size(paramCombinations,1)
    [logLike_all(i)] = fitFun(paramCombinations(i,:),parameterNames,psDataTr,1);
    if ~isServer, parfor_progress; end
end
if ~isServer, parfor_progress(0); end

% Finding the minimum of the log-likelihood values and the corresponding
% parameter combination
[logLike_min,idx] = min(logLike_all);
bestParamGs = paramCombinations(idx,:);
if ~isServer
    figure; plot(logLike_all,'k'); hold on; plot(idx,logLike_all(idx), 'r*')
end
% Printing elapsed time
cTemp = clock;
fprintf('Gridsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(cTemp,cStart)/86400,'dd HH:MM:SS'));
fprintf('Minimum negative log-likelihood: %.2f \n',logLike_min);

% Perform fminsearch to refine the best parameters
% Take the best parameters found in the gridsearch and use them as starting
% point in the fminsearch
fprintf('Performing fmincon\n');

% Options for fminsearch
opts = optimoptions('fmincon');
% opts = optimset('fminsearch');
if isServer, opts.Display = 'final'; else opts.Display = 'iter'; end
opts.TolFun = 1e-6;
opts.TolX = 1e-6;

LB = [cellfun(@min,gridVectors(1:end-1)),1];
UB = [cellfun(@max,gridVectors(1:end-1)),1];

% Creating anonymous function for input to fminsearchbnd
fun = @(param) fitFun(param,parameterNames,psDataTr,1);
[bestParamFm,logLikeFm] = fmincon(fun,bestParamGs,[],[],[],[],LB,UB,[],opts);
fprintf('fmincon elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(clock,cTemp)/86400,'dd HH:MM:SS'));
fprintf('Minimum negative log-likelihood: %.2f \n',logLikeFm);

% Evaluate fitted parameters to original hidden variables
[ll_est_train,mdEval_est_train] = fitFun(bestParamFm,parameterNames,psDataTr,1);
% Evalueate original parameters to original hidden variables
[ll_true_train,mdEval_true_train] = fitFun(trueParamVect,parameterNames,psDataTr,1);
% Evaluate fitted parameters generelized to new hidden variables
[ll_est_test,mdEval_est_test] = fitFun(bestParamFm,parameterNames,psDataTe,1);
% Evalueate original parameters generelized to new hidden variables
[ll_true_test,mdEval_true_test] = fitFun(trueParamVect,parameterNames,psDataTe,1);

% generating response and predicted distributions
funKs = @(x) ksdensity(x,'bandwidth',1);
    function out = wrapKs(x)
        [ff,fx] = funKs(x);
        out = {ff,fx};
    end
respDistr_train = varfun(@wrapKs,psDataTr,...
    'InputVariables',{'respA','respV'},'GroupingVariables',{'locV','locA'});
respDistr_train.f_respA = respDistr_train.Fun_respA(:,1);
respDistr_train.x_respA = respDistr_train.Fun_respA(:,2);
respDistr_train.f_respV = respDistr_train.Fun_respV(:,1);
respDistr_train.x_respV = respDistr_train.Fun_respV(:,2);
respDistr_train.Fun_respA = [];
respDistr_train.Fun_respV = [];
respDistr_train.GroupCount = [];
respDistr_train.Properties.RowNames = {};

respDistr_test = varfun(@wrapKs,psDataTe,...
    'InputVariables',{'respA','respV'},'GroupingVariables',{'locV','locA'});
respDistr_test.f_respA = respDistr_test.Fun_respA(:,1);
respDistr_test.x_respA = respDistr_test.Fun_respA(:,2);
respDistr_test.f_respV = respDistr_test.Fun_respV(:,1);
respDistr_test.x_respV = respDistr_test.Fun_respV(:,2);
respDistr_test.Fun_respA = [];
respDistr_test.Fun_respV = [];
respDistr_test.GroupCount = [];
respDistr_test.Properties.RowNames = {};

prDistr_est_train = cat(1,mdEval_est_train.conditions);
[prDistr_est_train.f_predV,prDistr_est_train.x_predV] = cellfun(funKs,{mdEval_est_train.sV_resp}','UniformOutput',false);
[prDistr_est_train.f_predA,prDistr_est_train.x_predA] = cellfun(funKs,{mdEval_est_train.sA_resp}','UniformOutput',false);

prDistr_true_train = cat(1,mdEval_true_train.conditions);
[prDistr_true_train.f_predV,prDistr_true_train.x_predV] = cellfun(funKs,{mdEval_true_train.sV_resp}','UniformOutput',false);
[prDistr_true_train.f_predA,prDistr_true_train.x_predA] = cellfun(funKs,{mdEval_true_train.sA_resp}','UniformOutput',false);

prDistr_est_test = cat(1,mdEval_est_test.conditions);
[prDistr_est_test.f_predV,prDistr_est_test.x_predV] = cellfun(funKs,{mdEval_est_test.sV_resp}','UniformOutput',false);
[prDistr_est_test.f_predA,prDistr_est_test.x_predA] = cellfun(funKs,{mdEval_est_test.sA_resp}','UniformOutput',false);

prDistr_true_test = cat(1,mdEval_true_test.conditions);
[prDistr_true_test.f_predV,prDistr_true_test.x_predV] = cellfun(funKs,{mdEval_true_test.sV_resp}','UniformOutput',false);
[prDistr_true_test.f_predA,prDistr_true_test.x_predA] = cellfun(funKs,{mdEval_true_test.sA_resp}','UniformOutput',false);

% Neatly organizing fake data
pseudoDataPerCond = varfun(@(x) {x},psDataTr,...
    'InputVariables',{'respA','respV'},'GroupingVariables',{'locV','locA'});
pseudoDataPerCond.Properties.RowNames = {};
pseudoDataPerCond.Properties.VariableNames{'Fun_respA'} = 'respA';
pseudoDataPerCond.Properties.VariableNames{'Fun_respV'} = 'respV';
pseudoDataPerCond.GroupCount = [];

% % Removing unnecessary fields from bciSimulations
% fieldsToRemove = {'sA_resp','sV_resp','sV_hat','sA_hat','s_hat_common',...
%     'sV_hat_indep','sA_hat_indep'};
% for i = 1:numel(fieldsToRemove)
%     if isfield(mdEval_fit_new,fieldsToRemove{i});
%         mdEval_fit_new = rmfield(mdEval_fit_new,fieldsToRemove{i});
%     end
%     if isfield(mdEval_orig_new,fieldsToRemove{i});
%         mdEval_orig_new = rmfield(mdEval_orig_new,fieldsToRemove{i});
%     end
% end

% Saving subject specific bci simulations
workDir = fullfile(DEC_2_setupdir('final','anal_eeg'),'tests');
savePath = fullfile(workDir,['bci_paramrecov_',datestr(now,'yymmddHHMMSS'),'.mat']);
save(savePath,'-v7.3');

end