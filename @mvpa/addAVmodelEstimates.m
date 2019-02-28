function mvparesObj = addAVmodelEstimates(mvparesObj,varargin)
% Method for computing AV  model estimaties
% 
% USAGE:
%   mvparesObj = addAVmodelEstimates(mvparesObj)
% INPUT:
%   Required: 
%       mvparesObj (object): mvpares object
%   'Name'-Value arguments:
%       factors (cell array of strings): 
%           
%       splitExamples (string): split examples according to certain 
%           conditions and compute AV model estimate separately for splits. 
%           Accepted values: 'prestim_alpha_power','prestim_beta_power',
%           'prestim_theta_power', 'none' (default)
%       
% OUTPUT:
%   mvparesObj (object): mvpares object with AV model estimates

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validSplits = {'prestim_alpha_power','prestim_beta_power',...
    'prestim_theta_power', 'none'};
validFactors = {'Disp','VR','Task'};
addRequired(p,'mvparesObj',@(x) isa(x,'mvpares') && x.isvalid);
addParameter(p,'factors',{'Disp','VR','Task'},@(x) all(ismember(x, ...
                                                  validFactors)));
addParameter(p,'splitExamples','none',@(x) any(validatestring(x, ...
                                                  validSplits)));
addParameter(p,'discretize',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));

parse(p,mvparesObj,varargin{:});
mvparesObj = p.Results.mvparesObj;
factors = p.Results.factors;
splitExamples = p.Results.splitExamples;
discretize = p.Results.discretize;

% Defining constant parameters for AV weight computation
cParam.smallDisp = mvpa.smallDisp;
cParam.relVlevels = mvpa.relVlevels;
cParam.ciAlpha = 0.32;

if isvector(mvparesObj.getGenTimePoints)
    genPredLabels = mvparesObj.getPredLabels;
else
    genPredLabels = mvparesObj.getPredLabels('genTime','tr_x_tr');
end

% Check whether there are AV trials and whether there are incongruent ones 
% as well
genCond = mvparesObj.info.gen_cond;
if isempty(regexp(genCond,'AV-(i|ci).*','once'))
    warning('mvpa:addAVmodelEstimates:unsupportedCondition',...
        ['Audio-visual and incongruent data are requred to be able to ',...
        'compute AV weights. The generalization condition does not ',...
        'satisfy these. Returning.']);
    return;
end

% Starting the timer and printing details
cStart = clock;
ds = datestr(now);
printFnTitle(80,'addAVmodelEstimates',ds)
fprintf('Estimating... \n');

genExamplesInfo = mvparesObj.getInfoGenExamples;
if isempty(genPredLabels) || isempty(genExamplesInfo)
    warning('mvpa:addAVmodelEstimates:requestedDataNotPresent',...
        ['Couldn''t extract predicted labels/info of generalization ',...
        'examples, AV model estimates can''t be computed, returning.']);
    return;
end

if ~strcmp(splitExamples,'none')
    switch splitExamples
        case 'prestim_alpha_power'
            varName = 'psAlphaPowOcc';
            splitLevels = {'psAlow','psAhigh'};
        case 'prestim_beta_power'
            varName = 'psBetaPowTemporoPar';
            splitLevels = {'psBlow','psBhigh'};
        case 'prestim_theta_power'
            varName = 'psThetaPowFrontal';
            splitLevels = {'psTlow','psThigh'};
    end     
    if ismember(varName,genExamplesInfo.Properties.VariableNames)
        m = nanmedian(genExamplesInfo.(varName));
        splits = ones(size(genExamplesInfo,1),1);
        splits(genExamplesInfo.(varName) > m) = 2;
    else
        warning('mvpa:addAVmodelEstimates:requestedDataNotPresent',...
            ['The dataset does not contain the ''%s'' variable, ',...
            'AV model estimates can''t be computed, returning.'],varName);
        return;
    end
else
    splits = ones(size(genExamplesInfo,1),1);
end
genExamplesInfo = genExamplesInfo(:,{'locA','locV','relV','task'});

% Preparing gridsearch parameters
if discretize
    responseLoc = unique(genExamplesInfo.locV(~isnan(genExamplesInfo.locV)));
    genPredLabels = arrayfun(@(x) mvpa.chooseClosest(responseLoc,x),genPredLabels);
end

% If there is no parallel pool running, open one.
currPool = gcp('nocreate');
if isempty(currPool)
    parpool('local');
end

nSplits = length(unique(splits));

for i = 1:nSplits
    actPredLabels = genPredLabels(splits == i,:,:,:);
    actExamplesInfo = genExamplesInfo(splits == i,:);
    % Computing betas
    betas = mvpa.fitparams(actPredLabels,actExamplesInfo,factors,cParam);
    
    % Setting the dataset object to writable if it is not
    if ~mvparesObj.writable, mvparesObj.setWritable(true); end
    if nSplits == 1
        fieldName = 'gen_AVmodelEstimates';
    else
        fieldName = ['gen_AVmodelEstimates','_',splitLevels{i}];
    end
    if ismember(fieldName,fieldnames(mvparesObj.data))
        warning('mvpa:addAVmodelEstimates:overwriteField',...
            ['The field %s already exists in the mvpa result dataset, '...
            'it will be overwritten.'],fieldName);
    end
    mvparesObj.data.(fieldName) = orderfields(betas);
    
end

mvparesObj.setWritable(false);
% Finishing timer and printing elapsed time
fprintf('Estimation elapsed time (days hours:minutes:seconds) %s \n\n',...
    datestr(etime(clock,cStart)/86400,'dd HH:MM:SS'));

end

 