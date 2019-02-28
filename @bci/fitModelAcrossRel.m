function varargout = fitModelAcrossRel(parameters,parameterNames,dataVA,responseLoc,model,varargin)

% Parsing input
p = inputParser;

validModels = {'bci','bciSel','bciMatch','fus','segA','segV','taskRel','null'};
validOutputs = {'logLike','full'};

addRequired(p,'parameters',@(x) validateattributes(x,{'numeric'},{'vector'}));
addRequired(p,'parameterNames',@(x) validateattributes(x,{'cell'},...
    {'numel',numel(parameters)}));
addRequired(p,'dataVA',@(x) validateattributes(x,{'table'},{'nonempty'}));
addRequired(p,'responseLoc',@(x) validateattributes(x,{'numeric'},{'vector'}));
addRequired(p,'model',@(x) any(validatestring(x,validModels)));
addOptional(p,'decisionFun',1,@(x) validateattributes(x,{'numeric'},...
    {'scalar','integer','positive','<=',3}));
addOptional(p,'output','logLike',@(x) any(validatestring(x,validOutputs)));

parse(p,parameters,parameterNames,dataVA,responseLoc,model,varargin{:});

parameters = p.Results.parameters;
parameterNames = p.Results.parameterNames;
dataVA = p.Results.dataVA;
responseLoc = p.Results.responseLoc;
model = p.Results.model;
decisionFun = p.Results.decisionFun;
output = p.Results.output;

% Making sure parameters and parameterNames are row vectors
if iscolumn(parameters), parameters = parameters'; end
if iscolumn(parameterNames), parameterNames = parameterNames'; end

% Visual reliability levels in the data
relLevels = unique(dataVA.relV);
nLevelsRel = numel(relLevels);

if ~ismember(model,{'segA','null'})
    sigVnames = parameterNames(~cellfun(@isempty,regexp(parameterNames,'sigV[0-9]*')));
    otherParamIdx = cellfun(@isempty,regexp(parameterNames,'sigV[0-9]*'));
    
    if nLevelsRel ~= numel(sigVnames)
        error('bci_fitmodel_across_rel:inconsistentInput',...
            ['The number of sigV parameters must match the number of reliablity ',...
            'levels specified in the data.']);
    end
end

negLogLike_sum = 0;
if strcmp(output,'full'), mdlEval = struct([]); end

% loop over number of reliability levels
for iRel = 1:nLevelsRel
    
    % Select data for each reliability level
    actDataVA = dataVA((dataVA.relV == relLevels(iRel)),:); 
    
    if ~ismember(model,{'segA','null'})
        % Select visual variance parameter for this reliability level
        actSigV = parameters(ismember(parameterNames,sigVnames{iRel}));
        % generate the parameter setting for one visual reliability level
        actParameters = [parameters(otherParamIdx),actSigV];
        actParamNames = [parameterNames(otherParamIdx),'sigV'];
    else
        actParameters = parameters;
        actParamNames = parameterNames;
    end
    
    switch output
        
        case 'logLike'
            
            switch model
                case 'bci'
                    negLogLike_act = bci.fitModel(actParameters,actParamNames,actDataVA,responseLoc,decisionFun);
                case 'bciSel'
                    negLogLike_act = bci.fitModel(actParameters,actParamNames,actDataVA,responseLoc,decisionFun);
                case 'bciMatch'
                    negLogLike_act = bci.fitModel(actParameters,actParamNames,actDataVA,responseLoc,decisionFun);
                case 'fus'
                    negLogLike_act = bci.fitModelFus(actParameters,actParamNames,actDataVA,responseLoc);
                case 'segA'
                    negLogLike_act = bci.fitModelSegA(actParameters,actParamNames,actDataVA,responseLoc);
                case 'segV'
                    negLogLike_act = bci.fitModelSegV(actParameters,actParamNames,actDataVA,responseLoc);
                case 'taskRel'
                    negLogLike_act = bci.fitModelTaskRel(actParameters,actParamNames,actDataVA,responseLoc);
                case 'null'
                    negLogLike_act = bci.fitModelNull(actDataVA,responseLoc);
            end
            
        case 'full'
            
            switch model
                case 'bci'
                    [negLogLike_act,mdlEval_act] = bci.fitModel(actParameters,actParamNames,...
                                                                actDataVA,responseLoc,decisionFun);
                case 'bciSel'
                    [negLogLike_act,mdlEval_act] = bci.fitModel(actParameters,actParamNames,...
                                                                actDataVA,responseLoc,decisionFun);
                case 'bciMatch'
                    [negLogLike_act,mdlEval_act] = bci.fitModel(actParameters,actParamNames,...
                                                            actDataVA,responseLoc,decisionFun);
                case 'fus'
                    [negLogLike_act,mdlEval_act] = bci.fitModelFus(actParameters,actParamNames,...
                                                                   actDataVA,responseLoc);
                case 'segA'
                    [negLogLike_act,mdlEval_act] = bci.fitModelSegA(actParameters,actParamNames,...
                                                                    actDataVA,responseLoc);
                case 'segV'
                    [negLogLike_act,mdlEval_act] = bci.fitModelSegV(actParameters,actParamNames,...
                                                                    actDataVA,responseLoc);
                case 'taskRel'
                    [negLogLike_act,mdlEval_act] = bci.fitModelTaskRel(actParameters,actParamNames,...
                                                                      actDataVA,responseLoc);
                case 'null'
                    error('bci:fitModelAcrossRel:invalidParameter',...
                          ['''Full'' output is not available with ' ...
                           'null model']);
            end
            
            temp = repmat({relLevels(iRel)},size(mdlEval_act));
            [mdlEval_act.relV] = temp{:};
            mdlEval = cat(2,mdlEval,mdlEval_act);
            
    end
    
    % sum logLikes over reliability levels
    negLogLike_sum = negLogLike_sum + negLogLike_act; 
    
end

varargout{1} = negLogLike_sum;
if strcmp(output,'full')
    % Computing Bayesian Information Criterion
    n = size(dataVA,1);
    m = numel(parameterNames);
    bic = (-negLogLike_sum)-(0.5*m*log(n));
    varargout{2} = bic;
    varargout{3} = mdlEval;
end

end


        
