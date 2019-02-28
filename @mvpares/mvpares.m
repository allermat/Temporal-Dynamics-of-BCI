classdef mvpares < handle
    
    properties (SetAccess = private)
        level
        plotTimeWin
        source
        state
        writable
    end
    
    properties (Access = {?mvpares,?mvpa})
        data
        info
    end
    
    methods
        
        function obj = mvpares(fileName)
            % Class for mvpa results
            % 
            % DETAILS: 
            %   This class implements all methods for accessing and
            %   plotting various pieces of data from the mvpa result
            %   datasets. For creating and manipulating mvpa result
            %   datasets see the mvpa class. 
            % USAGE: 
            %   obj = mvpares(fileName)
            % INPUT:
            %   fileName (string): name of the mvpa result dataset from
            %       which the object is constructed.
            % OUTPUT: 
            %   obj (object): mvpares object
            
            % Copyright(C) 2016, Mate Aller
            % allermat@gmail.com
            
            % parsing input
            p = inputParser;
            addRequired(p,'fileName',@exist);
            parse(p,fileName);
            fileName = p.Results.fileName;
            
            % Initializing variables
            mResData = matfile(fileName);
            mResInfo = mResData.info;
            obj.data = mResData;
            obj.info = mResInfo;
            if strcmp(mResInfo.subID,'group')
                obj.level = 'group';
                obj.state = 'trained_and_generalized';
            else
                obj.level = 'subj';
                if all(ismember({'gen_accuracies','gen_groupings',...
                        'gen_examples','gen_predlabels'},fieldnames(mResData)));
                    obj.state = 'trained_and_generalized';
                else
                    obj.state = 'trained';
                end
            end
            obj.plotTimeWin = [mean(mResInfo.tr_timePoints{1}),mean(mResInfo.tr_timePoints{end})];
            obj.source = mResData.Properties.Source;
            obj.writable = false;
            
        end
        
        % Method for accessing the info property
        out = getInfo(obj)
        
        % Method for accessing AV model estimates
        corrCoeffs = getAVmodelCorrelations(obj,var,varargin)
        
        % Method for accessing AV model estimates
        estimats = getAVmodelEstimates(obj,varargin)
        
        % Method for accessing AV model weights
        weigths = getAVmodelWeights(obj,varargin)
        
        % Method for accessing AV model weights statistics
        stats = getAVmodelWeightStats(obj,varargin)
        
        % Method for accessing the trained weights for features
        w = getFeatureWeights(obj)
        
        % Get the training sampling frequency
        fs = getFsample(obj)
        
        % Method for getting prediction performance estimates
        estimates = getGenPerfEstimates(obj,varargin);
        
        % Method for getting the generalization time points
        genTimePoints = getGenTimePoints(obj)
        
        % Method for accessing the table of various information about the 
        % generalization examples 
        genExamplesInfo = getInfoGenExamples(obj,varargin)
        
        % Method for getting predicted labels
        predLabels = getPredLabels(obj,varargin)
        
        % Method for getting true labels
        trueLabels = getTrueLabels(obj,varargin)
        
        % Method for getting the size of the generalization time points
        s = getSizeGenTime(obj)
        
        % Method for getting the size of the generalization time points
        s = getSizeTrTime(obj)
        
        % Method for getting stats for various variables
        stats = getStats(obj,varName,varargin)
        
        % Method for getting the training time points
        trTimePoints = getTrTimePoints(obj)
        
        % Method for finding the index of a particular training time point
        idx = indTrTimePoint(obj,timePoint)
        
        % Method for setting the time window for plotting
        setPlotTimeWin(obj,timeWin)
        
        % Method for setting the 'writable' property
        setWritable(obj,isWritable)
        
        % Method for plotting AV model correlations
        hFig = showAVmodelCorrelations(obj,var,varargin)
        
        % Method for plotting estimated AV model betas
        hFig = showAVmodelEstimates(obj,effect,varargin)
        
        % Method for plotting estimated AV model weigths
        hFig = showAVmodelWeights(obj,effect,varargin)
        
        % Method for plotting prediction performance estimates
        hFig = showGenPerf(obj,varargin)
        
        % Method for plotting predicted labels vs. true labels
        hFig = showPredLabels(obj,varargin)
        
        % Method to keep the functionality of the Matlab.io.matFile
        % object's who method
        out = who(obj)
        
        % Method to keep the functionality of the Matlab.io.matFile
        % object's whos method
        out = whos(obj)
        
    end
    
    methods (Access = {?mvpares,?mvpa})
        
        % Method for getting the number of cross-validation folds
        out = getNcvFolds(obj)
        
        % Method for determining if the training condition and the
        % generalization condition are identical.
        out = getTrGenCondMatch(obj)
        
    end
    
    methods (Static)
        
        % Method for finding clusters of true values in input vector
        clusters = findClusters(h)
        
        % Method for plotting distributions of predicted labels
        plotLabelDistribution(h,data,groups,xLabelStr,yLabelStr,titleStr,y_lim)
        
        % Method for plotting data of time x time matrix format
        [hFig,hAxes] = plotTimeByTimeMatrix(time,data,varargin)
        
        % Method for plotting data of timeseries format
        [hFig,hAxes] = plotTimeSeriesData(time,data,varargin)
        
        % Method for segmenting figure area for multiple plots
        [nRowSubplot,nColSubplot,subPlotIdx] = segmentFigArea(nRows,nCols,spacing)
        
    end
    
end