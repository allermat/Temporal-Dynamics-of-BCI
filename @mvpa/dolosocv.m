function [models,scParam,acc,cvMisc] = dolosocv(feats,labs,grouping,cvDetails)
% Performs a leave one session out crossvalidation

% Parsing input
p = inputParser;

reqFieldsCvDetails = {'doGridSearch','maxnSVs','nFeatures',...
    'nModelsPerTimePoint','nParam','params','scMethod','svmType'};

% Defining inputs.
addRequired(p,'feats');
addRequired(p,'labs');
addRequired(p,'grouping');
addRequired(p,'cvDetails',@(x) all(isfield(x,reqFieldsCvDetails)));

% Parsing inputs.
parse(p,feats,labs,grouping,cvDetails);

% Assigning inputs to variables. 
doGridSearch = p.Results.cvDetails.doGridSearch;
feats = p.Results.feats;
grouping = p.Results.grouping;
labs = p.Results.labs;
cvDetails = p.Results.cvDetails;
maxnSVs = cvDetails.maxnSVs;
nFeatures = cvDetails.nFeatures;
nModelsPerTimePoint = cvDetails.nModelsPerTimePoint;
nParam = cvDetails.nParam;
params = cvDetails.params;
scMethod = cvDetails.scMethod;
svmType = cvDetails.svmType;

% Initializing arrays for data collection (padded with NaNs).
if any(ismember(svmType,{'cc','nc'}))
    nClassesMdl = max([2,numel(unique(labs))]);
else
    nClassesMdl = 2;
end
models = NaN(nFeatures+nClassesMdl+1,maxnSVs,nModelsPerTimePoint);
scParam = NaN(nFeatures,2,nModelsPerTimePoint);
acc = NaN(3,nModelsPerTimePoint);

if doGridSearch
    if ~isfield(cvDetails,'nStepsFine');
        error('Missing variables ''nStepsFine'' in cvDetails! ');
    else
        nStepsFine = cvDetails.nStepsFine;
    end
    if nParam == 1
        gridsCoarse = NaN(size(params{1},2),nModelsPerTimePoint);
        gridsFine = NaN(nStepsFine,nModelsPerTimePoint);
    elseif nParam == 2
        gridsCoarse = NaN(size(params{1},2),size(params{2},2),nModelsPerTimePoint);
        gridsFine = NaN(nStepsFine,nStepsFine,nModelsPerTimePoint);
    end
    paramFine1 = NaN(nStepsFine,nModelsPerTimePoint);
    paramFine2 = NaN(nStepsFine,nModelsPerTimePoint);
    bestParCoarse = squeeze(NaN(nParam,nModelsPerTimePoint));
    bestParFine = squeeze(NaN(nParam,nModelsPerTimePoint));
end

for iModel = 1:nModelsPerTimePoint
    
    % Scaling features separately for each model
    if strcmp(scMethod,'Z-f')
        trFeats = feats(grouping ~= iModel,:);
    else
        [trFeats,scParam(:,:,iModel)] = mvpa.scalefeatures(feats(grouping ~= iModel,:),scMethod);
    end
    trLabs = labs(grouping ~= iModel);
    
    if doGridSearch
        [finalParam,nCVdetails] = mvpa.gridsearch(trLabs,trFeats,svmType,...
            params,nModelsPerTimePoint-1,nStepsFine);
        if nParam == 1
            gridsCoarse(:,iModel) = nCVdetails{1,1};
            gridsFine(:,iModel) = nCVdetails{2,1};
            paramFine1(:,iModel) = nCVdetails{2,2};
            bestParCoarse(iModel) = nCVdetails{1,4};
            bestParFine(iModel) = nCVdetails{2,4};
        elseif nParam == 2
            gridsCoarse(:,:,iModel) = nCVdetails{1,1};
            gridsFine(:,:,iModel) = nCVdetails{2,1};
            paramFine1(:,iModel) = nCVdetails{2,2};
            paramFine2(:,iModel) = nCVdetails{2,3};
            bestParCoarse(:,iModel) = nCVdetails{1,4};
            bestParFine(:,iModel) = nCVdetails{2,4};
        end
    else
        finalParam(1) = params{1};
        if nParam == 2, finalParam(2) = params{2}; end
    end
    
    if strcmp(svmType,'cc')
        args = sprintf('-s 0 -t 0 -c %g -q',2^finalParam(1));
    elseif strcmp(svmType,'nc')
        args = sprintf('-s 1 -t 0 -n %g -q',finalParam(1));
    elseif strcmp(svmType,'er')
        args = sprintf('-s 3 -t 0 -c %g -p %g -q',2^finalParam(1),2^finalParam(2));
    elseif strcmp(svmType,'nr')
        args = sprintf('-s 4 -t 0 -c %g -n %g -q',2^finalParam(1),finalParam(2));
    end
    
    % Training model
    mdl = svmtrain(trLabs,trFeats,args);
    % Predicting using the training set to get the training error
    [~,acc(:,iModel),~] = svmpredict(trLabs,trFeats,mdl);
    
    % Collecting data.
    mdlmat = mvpa.mdl2mat(mdl);
    models(1:size(mdlmat,1),1:size(mdlmat,2),iModel) = mdlmat;
    cvMisc = struct();
    if doGridSearch
        cvMisc.gridsCoarse = gridsCoarse;
        cvMisc.gridsFine = gridsFine;
        cvMisc.paramFine1 = paramFine1;
        if nParam == 2
            cvMisc.paramFine2 = paramFine2;
        end
        cvMisc.bestParCoarse = bestParCoarse;
        cvMisc.bestParFine = bestParFine;
    end
    
end

end