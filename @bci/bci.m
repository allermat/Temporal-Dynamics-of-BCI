classdef bci
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        % Number of internally generated samples for fitting
        nIntSamples = 10000;
    end
    
    methods (Static)
        
        % Find the N best performing parameter combinations weighted by their Eucledian distance
        idx = chooseNbestWeighted(param,perf,n,thres)
        
        % Converts predicted labels to the format accepatble for the modelfitting
        bciData = prepNeuralData(info,actSample,varargin)
        
        % Fit full Bayesian Causal Inference model
        [logLike,det] = fitModel(parameters,parameterNames,dataVA,responseLoc,decisionFun)
        
        % Wrapper function for fitting across multiple reliability levels
        varargout = fitModelAcrossRel(parameters,parameterNames,dataVA,responseLoc,model,varargin)
        
        % Simulates responses from the forced fusion model for a given parameter set
        [logLike,det] = fitModelFus(parameters,parameterNames,dataVA,responseLoc)
        
        % Simulates responses from the segregation auditory model for a given parameter set
        [logLike,det] = fitModelSegA(parameters,parameterNames,dataVA,responseLoc)
        
        % Simulates responses from the segregation visual model for a given parameter set
        [logLike,det] = fitModelSegV(parameters,parameterNames,dataVA,responseLoc)
        
        % Simulates responses from the task relevance model for a given parameter set
        [logLike,det] = fitModelTaskRel(parameters,parameterNames,dataVA,responseLoc)
        
        % Simulates responses from a null model
        logLike = fitModelNull(parameters,parameterNames,dataVA,responseLoc)
        
        % Metohd for generating fake responses according to different models
        resp = generateFakeResponses(cfg)
        
        % Method for running parameter recovery for a specified model and parameters
        paramRecovery(I)
        
        % Method for plotting detailed results from parameter recovery
        hFig = plotRecovDtld(respDistr,predDistr1,predDistr2,varargin)
        
        % Method for plotting true vs predictied continuous distributions
        hFig = plotRespVsPredDistr(respDistr,predDistri)
        
    end
    
end

