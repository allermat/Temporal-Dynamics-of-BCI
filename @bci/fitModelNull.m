function negLogLike = fitModelNull(dataVA,responseLoc)
% Generates responses from a null model for a given parameter set
% 
% DETAILS: 
%   A uniform response distribution at chance level is
%   created for every condition
% USAGE:
%   [logLike,all] = bci.fitModelNull(dataVA,responseLoc)
% INPUT:
%   dataVA (table): behavioural response data with which the 
%       predicted responses of the BCI model will be compared. Given in
%       visual degrees. 
%       Rows = trials, Variables: locA, locV, respA, respV
%   responseLoc is the locations corresponding to the response bottons, e.g.
%       degrees. If a single value, continuous response with a
%       kernel width is used.
% OUTPUT:
%   negLogLike:
%   mdlEval: 

% This code was contributed by Ulrik Beierholm (Apr 2016), with some code
% left over from Wei Ji Ma and Konrad Koerding at some time 2007-2009
% Further expanded and adjusted by Uta Noppeney & Mate Aller

% Fix the initialising random state - WHY ARE WE USING A FIXED SEED?
% rand('state',55); 
% randn('state',55);
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

% discrete responses?? i.e. response loc vector e.g. [1 2 3 4 5]
% or continuous responseLoc = 1
if length(responseLoc) == 1
    error('bci:fitmodelNull:notImplemented',...
          ['Null model is not yet implemented for the continuous ' ...
           'case']);
end

% Initialze variable to collect log-likelihood values for each condition
logLikeCond = NaN(nCond,1);

% Simulate responses for each condition
for indCond = 1:nCond
    
    % A, V responses given by participant for particular A,V location
    % combination
    dataV = dataConditions{indCond}(:,1)'; dataV = dataV(~isnan(dataV));
    dataA = dataConditions{indCond}(:,2)'; dataA = dataA(~isnan(dataA));
    
    %  Compute loglike for A and V responses

    % Uniform distribution at chance level
    freq_predV = (1/numel(responseLoc))*ones(1,numel(responseLoc));
    freq_predA = freq_predV;
    
    %calculate absolute frequencies of actual responses
    freq_dataV = hist(dataV, responseLoc);
    freq_dataA = hist(dataA, responseLoc);
    
    %calculate log-likelihood
    logLikeA = sum(freq_dataA .* log(freq_predA)) ;
    logLikeV = sum(freq_dataV .* log(freq_predV)) ;

    
    %---------------------------------------------------------------------
    % sum loglike across different task responses
    if all(isnan(logLikeA)) && all(isnan(logLikeV))
        logLikeCond(indCond) = NaN;
    else
        % If all(isnan(logLikeA)) || all(isnan(logLikeV)), then the sum 
        % is going to be 0 whis means, the ligelihood is 1, which is not
        % true. 
        logLikeCond(indCond) = nansum([nansum(logLikeA) nansum(logLikeV)]);
    end
    
end  % end of loop over conditions

%sum over conditions and turn into negative log likelihood
negLogLike = -sum(logLikeCond); % neg log like for each reliability level
    
end