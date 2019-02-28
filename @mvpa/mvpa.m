classdef mvpa

    properties (Constant)
        
        % Small AV disparity in degrees
        smallDisp = 6.67;
        % Visual reliability levels
        relVlevels = [2,12];
        
    end
    
    methods (Static)
        
        % Method for computing correleations of AV weights
        mvparesObj = addAVmodelCorrelations(mvparesObj,varargin)
        
        % Method for computing AV model estimaties
        mvparesObj = addAVmodelEstimates(mvparesObj,varargin)
        
        % Method for computing generalization performance estimates
        mvparesObj = addGenPerfEstimates(mvparesObj)
        
        % Method for computing various statistics
        mvparesObj = addStats(mvparesObj,var,varargin)
        
        % Method for appending data to an existing MVPA dataset
        appenddata(I)
        
        % Averages mvpares objects
        mvparesObj = averageMvpaRes(I)
        
        % Method for choosing array entry closest to number x
        c = chooseClosest(vals,x)
            
        % Method for merging mvpares objects
        mvparesObj = mergeMvpaRes(I)
        
        % Prepares an MVPA dataset and saves it on hard drive
        prepdata(I)
        
        % Method for plotting EEG topographies from the MVPA dataset
        showTopographies(filePath)
        
        % Method for generalizing a support vector machine model on
        % timeseries data
        [mvparesObj,filePath] = svmGeneralize(I)
        
        % Method for training a support vector machine model on timeseries 
        % data
        [mvparesObj,filePath] = svmTrain(I)
        
    end
    
    methods (Static) %, Access = private) 
        
        % Method for computing confidence interval of mean
        [ciLow,ciHigh] = confmean(A,dim,alpha)
        
        % Method for groups the examples for cross-validation
        groupings = assigngroupings(cvScheme,trExamples_info,misc)
        
        % Method to compute statistics for AV model correlations.
        stats = compAVmodelCorrStats(indivData,time,genTime,varargin)
        
        % Method to compute statistics for AV weights data.
        stats = compAVmodelWeightStats(indivData,time,varargin)
        
        % Method to compute statistics for decoding performance
        stats = compGenPerfStats(indivData,time,varargin)
        
        % Method to compute various performance estimates across
        % training and generalization samples.
        varargout = compPerfEstimates(trueLab,predLab,mode)
        
        % Method for performing a full K-fold crossvalidation
        [models,scParam,acc,cvMisc] = dokfoldcv(feats,labs,groupings,cvDetails)
        
        % Method for performing a leave one session out crossvalidation
        [models,scParam,acc,cvMisc] = dolosocv(feats,labs,grouping,cvDetails)
        
        % Extracts the necessary data for MVPA for a given subject.
        [feat,info,misc] = extractdata(dataDir,condDef,idStr,tr_method)
        
        % Extracts the required set of features
        feat = extractfeatures(dataM,isExample,trMethod,trLabel,timePoints,varargin)
        
        % Method for fitting design parameters to predicted labels
        betas = fitparams(genPredLabels,genExamples,smallDiscr,relV,ciAlpha)
        
        % Method for fitting design parameters to predicted labels
        betas = fitparams_pseu(genPredLabels,genExamples,smallDiscr,relV,ciAlpha)
        
        % Method for generating corrected AV data
        [info,infoCorr,featCorr] = generateCorrectedAvData(info,feat,condDef,corrMode,dataType)
        
        % Method for generating pseudo-AV data
        [infoPseu,featPseu,groupingsPseu] = generatePseudoAvData(info,feat,groupings,inputCond,condDef)
        
        % Performs grid search with two passes on the given parameters with nested cross-validation.
        [bestparam,details] = gridsearch(TR_labs,TR_feats,svm_type,paramCoarse,nFold,nStepsFine)
        
        % Returns the sample closest to some time point in the specified time point set.
        res = indsamplecustom(t,nsample,fs,starttime)
        
        % Converts the libsvm's svmtrain output model converted to a matrix back to sturcture form.
        mdl = mat2mdl(M)
        
        % Converts the libsvm's svmtrain output model to a matrix.
        M = mdl2mat(mdl)
                
        % Scales the features according to the specified method.
        [scaledFeats,varargout] = scalefeatures(feats,method,varargin)
        
        % Selects examples matching the given input condition
        isExample = selectexamples(inputConds,condDef,dataInfo)
        
        % Selects the label for the given input label and condition
        labelOut = selectlabel(inputLabel,inputConds)
        
    end
    
end