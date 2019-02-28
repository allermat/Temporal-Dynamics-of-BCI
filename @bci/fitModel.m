function [negLogLike,mdlEval] = fitModel(parameters,parameterNames,dataVA,responseLoc,decisionFun)
% Simulates responses from the BCI model for a given parameter set
% 
% DETAILS: 
%   Creating 10000 (e.g.) internal samples and creating a distribution of 
%   potential responses based on these. These distributions are then 
%   compared with the observed distributions during estimation to compute 
%   loglike for this parameter setting.
% USAGE:
%   [logLike,all] = bci.fitModel(parameters,parameterNames,dataVA,responseLoc,readouttype)
% INPUT:
%   parameters (numeric vector): BCI model parameters. I would rather have 
%       used a table instead of using parameters and parameterNames as 
%       separate inputs, but this input has to be compatible with
%       fminsearch. 
%   parameterNames (cell array of strings): BCI model parameter names.
%       Possible values: 'p_common','muP','sigP','kW','sigA','sigV'
%   dataVA (table): behavioural response data with which the 
%       predicted responses of the BCI model will be compared. Given in
%       visual degrees. 
%       Rows = trials, Variables: locA, locV, respA, respV
%   responseLoc is the locations corresponding to the response bottons, e.g.
%       degrees. If a single value, continuous response with a
%       kernel width is used.
%   decisionFun = 1: model averaging, 2: model selection, 3: probability
%       matching
% OUTPUT:
%   negLogLike:
%   mdlEval: 

% This code was contributed by Ulrik Beierholm (Apr 2016), with some code
% left over from Wei Ji Ma and Konrad Koerding at some time 2007-2009
% Further expanded and adjusted by Uta Noppeney & Mate Aller

% Fix the initialising random state
rng(55);

nIntSamples = bci.nIntSamples;  %% minimum required for model fitting

% Checking dataVA
reqDataVars = {'locV','locA','respV','respA'};
if any(~ismember(reqDataVars,dataVA.Properties.VariableNames))
    error('bci_fitmodel:missingDataVariable',...
        'Missing variable from dataVA in bci_fitmodel.')
end

% Ensure that only data are included where both A and V signals are used
dataVA = dataVA(~isnan(dataVA.locV) & ~isnan(dataVA.locA),:);

% Check responseLoc
if iscolumn(responseLoc)
    responseLoc = responseLoc';
end

% Check parameters
if iscolumn(parameters)
    parameters = parameters';
end

% To save computation time, check for what stimulus conditions were
% presented, then simulate each one

% Ensure, that the order of variables in dataVA is correct, also convert
% from table to matrix
dataVA = table2array(dataVA(:,{'locV','locA','respV','respA'}));
% Since condition labels are floating point numbers, we have to use a
% tolerance value for comparison operations. The unique function is not
% capable of this. The consolidator function from Matlab Central takes care
% of this operation, however consider swithcing to uniqetol if you have
% MATLAB R2015a or higher. 
tol = eps('single');
conditions = consolidator(dataVA(:,[1,2]),[],[],tol);
nCond = size(conditions,1);
dataConditions = cell(nCond,1);
for i = 1:nCond
    % Finding the examples corresponding to a condition. Since condition
    % labels are floating point numbers we have to use a tolerance value
    % for comparison operations 
    match = all(abs(dataVA(:,[1,2])-repmat(conditions(i,:),size(dataVA,1),1)) < tol,2);
    dataConditions{i} = dataVA(match,[3,4]);
end

% Parameters: [p_common(i),xP(m),sigP(m),kernelWidth(n),sigA(p),sigV(k)]
p_common = parameters(ismember(parameterNames,'p_common')); % Set to 1 for full integration, to 0 for full segregation
sigP = parameters(ismember(parameterNames,'sigP'));
sigA = parameters(ismember(parameterNames,'sigA'));
sigV = parameters(ismember(parameterNames,'sigV'));

% Setting default for muP if it is not specified
if any(ismember(parameterNames,'muP'))
    xP = parameters(ismember(parameterNames,'muP'));
else
    xP = 0;
end

% discrete responses?? i.e. response loc vector e.g. [1 2 3 4 5]
% or continuous responseLoc = 1
if length(responseLoc) == 1
    kernelWidth = parameters(ismember(parameterNames,'kW'));  %% kernel width as a parameter that is fitted to the data
else
    kernelWidth = 1;  % set to 1 for discrete responses
end

% Throw an error if there is a missing parameter
if any(cellfun(@isempty,{p_common,xP,sigP,sigA,sigV,kernelWidth}))
    error('bci_fitmodel:missingParameter','Missing parameter in bci_fitmodel.')
end

% Variances of A and V and prior
varV = sigV^2;
varA = sigA^2;
varP = sigP^2;

% Variances of estimates given common or independent causes
varVA_hat = 1/(1/varV + 1/varA + 1/varP);
varV_hat = 1/(1/varV + 1/varP);
varA_hat = 1/(1/varA + 1/varP);

% Variances used in computing probability of common or independent causes
var_common = varV * varA + varV * varP + varA * varP;
varV_indep = varV + varP;
varA_indep = varA + varP;

% Initialze variable to collect log-likelihood values for each condition
logLikeCond = NaN(nCond,1);

% Simulate responses for each condition
for indCond = 1:nCond
    
    sV = conditions(indCond,1);%sourcevec(sVind);
    sA = conditions(indCond,2);%sourcevec(sAind);
    
    % Generation of fake data
    xV = sV + sigV * randn(nIntSamples,1);
    xA = sA + sigA * randn(nIntSamples,1);
    
    % Estimates given common or independent causes
    s_hat_common = (xV/varV + xA/varA + xP/varP) * varVA_hat;
    sV_hat_indep = (xV/varV + xP/varP) * varV_hat;
    sA_hat_indep = (xA/varA + xP/varP) * varA_hat;
    
    % Probability of common or independent causes
    quad_common = (xV-xA).^2 * varP + (xV-xP).^2 * varA + (xA-xP).^2 * varV;
    quadV_indep = (xV-xP).^2;
    quadA_indep = (xA-xP).^2;
    
    % Likelihood of observations (xV,xA) given C, for C=1 and C=2
    likelihood_common = exp(-quad_common/(2*var_common))/(2*pi*sqrt(var_common));
    likelihoodV_indep = exp(-quadV_indep/(2*varV_indep))/sqrt(2*pi*varV_indep);
    likelihoodA_indep = exp(-quadA_indep/(2*varA_indep))/sqrt(2*pi*varA_indep);
    likelihood_indep =  likelihoodV_indep .* likelihoodA_indep;
    
    % Posterior probability of C given observations (xV,xA)
    post_common = likelihood_common * p_common;
    post_indep = likelihood_indep * (1-p_common);
    pC = post_common./(post_common + post_indep);
    
    % ----------------------------------------------------------------------
    % Generate spatial location responses and compute loglike
    switch decisionFun
        
        case 1
            % Mean of posterior - Model averaging
            % Overall estimates: weighted averages
            sV_hat = pC .* s_hat_common + (1-pC) .* sV_hat_indep;
            sA_hat = pC .* s_hat_common + (1-pC) .* sA_hat_indep;
            
        case 2
            % Model selection instead of model averaging
            % averaging
            sV_hat = (pC>0.5).*s_hat_common + (pC<=0.5).*sV_hat_indep;
            sA_hat = (pC>0.5).*s_hat_common + (pC<=0.5).*sA_hat_indep;
            
        case 3
            % Probability matching
            thresh = rand(nIntSamples,1); % to enable probability matching
            sV_hat = (pC>thresh).*s_hat_common + (pC<=thresh).*sV_hat_indep;
            sA_hat = (pC>thresh).*s_hat_common + (pC<=thresh).*sA_hat_indep;
            
    end
    
    % compute predicted responses for discrete and continuous case
    if length(responseLoc) > 1
        % discrete responses
        lengthRespLoc=length(responseLoc);
        
        %find the response location closest to sV_hat and
        %sA_hat, ie with minimum deviation
        [~,tV] = min(abs(repmat(sV_hat,1,lengthRespLoc) - repmat(responseLoc,nIntSamples,1)),[],2);
        [~,tA] = min(abs(repmat(sA_hat,1,lengthRespLoc) - repmat(responseLoc,nIntSamples,1)),[],2);
        sV_pred_resp = responseLoc(tV);
        sA_pred_resp = responseLoc(tA);
    else
        %continuous responses used, no discretisation
        sV_pred_resp = sV_hat;
        sA_pred_resp = sA_hat;
    end
    
    % A, V responses given by participant for particular A,V location
    % combination
    dataV = dataConditions{indCond}(:,1)'; dataV = dataV(~isnan(dataV));
    dataA = dataConditions{indCond}(:,2)'; dataA = dataA(~isnan(dataA));
    
    %  Compute loglike for A and V responses
    if length(responseLoc) > 1
        %discrete case
        %calculate frequencies of predictions, at least 0.00001 (1 out of 100,000)
        freq_predV = max(0.00001,hist(sV_pred_resp, responseLoc)/nIntSamples);
        freq_predA = max(0.00001,hist(sA_pred_resp, responseLoc)/nIntSamples);
        
        %calculate absolute frequencies of actual responses
        freq_dataV = hist(dataV, responseLoc);
        freq_dataA = hist(dataA, responseLoc);
        
        %calculate log-likelihood
        logLikeA = sum(freq_dataA .* log(freq_predA)) ;
        logLikeV = sum(freq_dataV .* log(freq_predV)) ;
        
    else
        %continuous case
        %gaussian kernel distribution for each condition
        %for each condition average over gaussian likelihoods
        
        logLikeA = log(mean(normpdf(repmat(dataA,nIntSamples,1) ,repmat(sA_pred_resp,1,length(dataA)), kernelWidth)));
        logLikeV = log(mean(normpdf(repmat(dataV,nIntSamples,1) ,repmat(sV_pred_resp,1,length(dataV)), kernelWidth)));
        
    end
    
    %---------------------------------------------------------------------
    % sum loglike across different task responses
    if all(isnan(logLikeA)) && all(isnan(logLikeV))
        % In this case the sum would be 0 which would mean that the
        % likelihood would be 1, which is not, so we prevent that
        % from happpening. 
        logLikeCond(indCond) = NaN;
    else
        logLikeCond(indCond) = nansum([nansum(logLikeA) nansum(logLikeV)]);
    end
    %---------------------------------------------------------------------
    % for a particular parameter setting make plots and save biasses etc.
    if nargout > 1
        
        %store values
        mdlEval(indCond).sA_resp = sA_pred_resp;
        mdlEval(indCond).sV_resp = sV_pred_resp;
        mdlEval(indCond).sV_hat = sV_hat;
        mdlEval(indCond).sA_hat = sA_hat;
        mdlEval(indCond).s_hat_common = s_hat_common;
        mdlEval(indCond).sV_hat_indep = sV_hat_indep;
        mdlEval(indCond).sA_hat_indep = sA_hat_indep;
        mdlEval(indCond).conditions = array2table(conditions(indCond,:),'VariableNames',{'locV','locA'});
        mdlEval(indCond).logLikeA = logLikeA;
        mdlEval(indCond).logLikeV = logLikeV;
        mdlEval(indCond).parameters = array2table(parameters,'VariableNames',parameterNames);
        
        if length(responseLoc) > 1
            mdlEval(indCond).freq_predV = freq_predV;
            mdlEval(indCond).freq_predA = freq_predA;
        end
    end
    
end  % end of loop over conditions

% sum over conditions and turn into negative log likelihood
negLogLike = -sum(logLikeCond);
    
end
