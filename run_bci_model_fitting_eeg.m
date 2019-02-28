function run_bci_model_fitting_eeg(subjListStr,varargin)

% Parsing input
p = inputParser;

addRequired(p,'subjListStr',@(x)validateattributes(x,{'char'},...
    {'nonempty'}));
addParameter(p,'discretize',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));
addParameter(p,'removeOutliers',false,@(x)validateattributes(x,{'logical'},...
    {'scalar'}));

parse(p,subjListStr,varargin{:});

subjListStr = p.Results.subjListStr;
discretize = p.Results.discretize;
removeOutliers = p.Results.removeOutliers;

% Basic parameters
subjList = regexp(subjListStr,'_','split');

expStage = 'final';
trMethod = 'sample-wise-sm-avg';

% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
    progrMonitor = false;
    addpath(genpath(fullfile(DEC_2_setupdir(expStage,'utils'),'Utility')));
else
    isServer = false;
    progrMonitor = true;
    dbstop if error;
end

% Opening parallel pool.
% if there is no parallel pool running, open one.
currPool = gcp('nocreate');
if isempty(currPool)
    if isServer
        parpool('local',16);
    else
        parpool('local');
    end
end

% function for generating response and predicted distributions
funKs = @(x) ksdensity(x,'bandwidth',1);
    function out = funKsCatOutput(x)
        [ff,fx] = funKs(x);
        out = {ff,fx};
    end
wrapKs = @(x) funKsCatOutput(x);

% % Function for choosing array entry closest to number x
%     function c = chooseClosest(vals,x)
%     [~,idx] = min(abs(vals-x));
%     c = vals(idx);
%     end

%% Loading data

fileMatchStr = 'er_tr-AxorV_stim_-100-5-1000_gen-AV-ci-av_stim_t_[0-9]+.mat';
% Printing header for console output
ds = datestr(now);
printFnTitle(80,'run_bci_model_fitting_eeg',ds)

for iSubj = 1:numel(subjList)
    subID = subjList{iSubj};
    
    % Starting the timer and printing details
    cStart = clock;
    fprintf('Processing subject %s... \n',subID);
    
    workDir = fullfile(DEC_2_setupdir(expStage,'anal_eeg_sub_mvpa',subID),trMethod);
    saveDf = cd(workDir);
    fileList = dir('*.mat');
    fileList = {fileList.name}';
    matchID = ~cellfun(@isempty,regexp(fileList,fileMatchStr));
    if sum(matchID) == 0
        warning('No file, skipping subject %s! ',subID);
        cd(saveDf);
        continue;
    elseif sum(matchID) > 1
        warning('Multiple files, skipping subject %s! ',subID);
        cd(saveDf);
        continue;
    else
        fileName = fileList(matchID);
        fileName = fileName{:};
        filePath = fullfile(workDir,fileName);
        cd(saveDf);
    end
    mvparesObj = mvpares(filePath);
    
    % Preparing eeg data
    predLabelInfo = mvparesObj.getInfoGenExamples;
    predLabel = mvparesObj.getPredLabels('genTime','tr');
    
    % Removing outliers. This is particularly important, when using
    % single trial data as the predicted labels can be noisy. I'm
    % using a standard cutoff, > 5 STD.
    if removeOutliers
        cutoff = 5;
        mu = mean(predLabel(:));
        s = std(predLabel(:));
        predLabel_Z = (predLabel-mu)/s;
        isOutl = any(abs(predLabel_Z) > cutoff, 2);
        nOutl = sum(isOutl); %#ok
        predLabel = predLabel(~isOutl,:);
        predLabelInfo = predLabelInfo(~isOutl,:);
    end
    
    % Preparing gridsearch parameters
    if discretize
        responseLoc = unique(predLabelInfo.locV(~isnan(predLabelInfo.locV)));
        predLabel = arrayfun(@(x) mvpa.chooseClosest(responseLoc,x),predLabel);
    else
        % Response locations, 1 if continuous response
        responseLoc = 1;
    end
    
    % Preparing the parameter space for gridsearch
    % grid resolution
    n = 5;
    p_common_full = linspace(0,1,n);
    sigV1 = linspace(0.1,30,n);
    sigV2 = linspace(0.1,30,n);
    sigA = linspace(0.1,30,n);
    sigP = linspace(0.1,30,n);
    kW = 1;
    
    % Generating parameter combinations for full model
    parameterNamesBCI = {'p_common','sigP','sigA','sigV1','sigV2','kW'};
    gridVectorsBCI = {p_common_full,sigP,sigA,sigV1,sigV2,kW};
    nParameters = numel(gridVectorsBCI);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectorsBCI{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinationsBCI = cat(2,coords{:});
    
    % Generating parameter combinations for fusion model
    parameterNamesFus = {'sigP','sigA','sigV1','sigV2','kW'};
    gridVectorsFus = {sigP,sigA,sigV1,sigV2,kW};
    nParameters = numel(gridVectorsFus);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectorsFus{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinationsFus = cat(2,coords{:});
    
    % Task relevance model
    parameterNamesTaskRel = {'sigP','sigA','sigV1','sigV2','kW'};
    paramCombinationsTaskRel = paramCombinationsFus;
    
    % Segregation visual model
    parameterNamesSegV = {'sigP','sigV1','sigV2','kW'};
    gridVectorsSegV = {sigP,sigV1,sigV2,kW};
    nParameters = numel(gridVectorsSegV);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectorsSegV{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinationsSegV = cat(2,coords{:});
    
    % Segregation visual model
    parameterNamesSegA = {'sigP','sigA','kW'};
    gridVectorsSegA = {sigP,sigA,kW};
    nParameters = numel(gridVectorsSegA);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectorsSegA{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinationsSegA = cat(2,coords{:});
    
    % Selecting time points of intersts
    timePoints = -0.1:0.005:0.7;
    fun = @(x) mvparesObj.indTrTimePoint(x);
    timePointIdx = arrayfun(fun,timePoints);
    predLabel = predLabel(:,timePointIdx);
    nSample = numel(timePoints);
    
    % Models to be fitted
    models = {'bci','fus','segA','segV','taskRel'};
    nModel = numel(models);
    
    % How to deal with the neural data selection from the cross-validation
    % folds. 
    neurDataSelection = 'none';
    
    % Pre-allocating arrays for data collection
    [negLogLike,negLogLikeAllGs,negLogLikeAllFm,bic,mdlEval,bestParamAllGs,bestParamFm,...
        bestParamAllFm,predDistr,R2] = deal(cell(nSample,nModel));
    respDistr = cell(nSample,1);
    paramCombinationsGs = {paramCombinationsBCI,paramCombinationsFus,...
        paramCombinationsSegA,paramCombinationsSegV,paramCombinationsTaskRel};
    parameterNames = {parameterNamesBCI,parameterNamesFus,...
        parameterNamesSegA,parameterNamesSegV,parameterNamesTaskRel};
    % Initializing progress monitor
    if progrMonitor, parfor_progress(nSample); end
    parfor iSample = 1:nSample
        
        bciData = bci.prepNeuralData(predLabelInfo,predLabel(:,iSample),neurDataSelection);
        
        for iModel = 1:nModel
            
            % Choosing model to fit
            actModel = models{iModel};
            
            % Decision function model averaging (1) vs. model selection (2)
            decisionFun = 1;
            
            switch actModel
                case 'bci'
                    paramCombinations = paramCombinationsBCI;
                    paramNames = parameterNamesBCI;
                case 'fus'
                    paramCombinations = paramCombinationsFus;
                    paramNames = parameterNamesFus;
                case 'segA'
                    paramCombinations = paramCombinationsSegA;
                    paramNames = parameterNamesSegA;
                case 'segV'
                    paramCombinations = paramCombinationsSegV;
                    paramNames = parameterNamesSegV;
                case 'taskRel'
                    paramCombinations = paramCombinationsTaskRel;
                    paramNames = parameterNamesTaskRel;
            end
            
            % Pre-allocating array for data collection
            logLike_temp = NaN(size(paramCombinations,1),1);
                        
            for i = 1:size(paramCombinations,1)
                [logLike_temp(i)] = bci.fitModelAcrossRel(paramCombinations(i,:),paramNames,bciData,responseLoc,actModel,decisionFun);
                
            end
            
            % Choose the best parameter combination(i.e. lowest logLike)
            % and other 10 best performing parameter combinations weighted
            % by their Eucledian distance, so they cover a broad range of
            % plausible parameter combinations
            thres = 1; %min([size(paramCombinations,1),100]);
            n = 1; %min([thres,10]);
            idx = bci.chooseNbestWeighted(paramCombinations,logLike_temp,n,thres);
            negLogLike_min = logLike_temp(idx(1));
            startParFm = paramCombinations(idx,:);
            
            bestParamAllGs{iSample,iModel} = startParFm;
            % if ~isServer
            %     figure();
            %     plot(logLike_temp,'k'); hold on;
            %     plot(idx,logLike_temp(idx),'r*');
            % end
            % Saveing all log-likelihoods and parameter combinations
            negLogLikeAllGs{iSample,iModel} = logLike_temp;
%             
%             % Printing elapsed time
%             cTemp = clock;
%             fprintf('Gridsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
%                 datestr(etime(cTemp,cStart)/86400,'dd HH:MM:SS'));
%             fprintf('Minimum negative log-likelihood: %.2f \n',logLike_min);
            
            %% Perform fminsearch to refine the best parameters
            % Take the best parameters found in the gridsearch and use them as starting
            % Options for fminsearch
            opts = optimset('fminsearch');
            if isServer, opts.Display = 'none'; else opts.Display = 'none'; end
            opts.MaxFunEvals = 4000;
            opts.TolFun = 1.e-6;
            opts.TolX = 1.e-6;
            opts.MaxIter = 1000;
%             opts = optimoptions('fmincon');
%             if isServer, opts.Display = 'none'; else opts.Display = 'final'; end
%             opts.TolFun = 1e-6;
%             opts.TolX = 1e-6;
            
            % set upper and lower bounds on parameters
            switch actModel
                case 'bci'
                    % 'p_common','sigP','sigA','sigV1','sigV2','kW'
                    LB = [cellfun(@min,gridVectorsBCI(1:end-1)),1];
                    UB = [cellfun(@max,gridVectorsBCI(1:end-1)),1];
                case 'fus'
                    % 'sigP','sigA','sigV1','sigV2','kW'
                    LB = [cellfun(@min,gridVectorsFus(1:end-1)),1];
                    UB = [cellfun(@max,gridVectorsFus(1:end-1)),1];
                case 'taskRel'
                    % 'sigP','sigA','sigV1','sigV2','kW'
                    LB = [cellfun(@min,gridVectorsFus(1:end-1)),1];
                    UB = [cellfun(@max,gridVectorsFus(1:end-1)),1];
                case 'segA'
                    % 'sigP','sigA','kW'
                    LB = [cellfun(@min,gridVectorsSegA(1:end-1)),1];
                    UB = [cellfun(@max,gridVectorsSegA(1:end-1)),1];
                case 'segV'
                    % 'sigP','sigV1','sigV2','kW'
                    LB = [cellfun(@min,gridVectorsSegV(1:end-1)),1];
                    UB = [cellfun(@max,gridVectorsSegV(1:end-1)),1];
            end
            
            bestParamTemp = NaN(size(startParFm));
            fvalTemp = NaN(size(startParFm,1),1);
            % Creating anonymous function for input to fminsearchbnd
            fun = @(param) bci.fitModelAcrossRel(param,paramNames,bciData,responseLoc,actModel,decisionFun);
            for i = 1:size(startParFm,1)
%                 [bestParamTemp(i,:),fvalTemp(i)] = ...
%                     fmincon(fun,startParFm(i,:),[],[],[],[],LB,UB,[],opts);
                [bestParamTemp(i,:),fvalTemp(i)] = ...
                    fminsearchbnd(fun,startParFm(i,:),LB,UB,opts);
            end
            
            [~,idx] = min(fvalTemp);
            bestParamFm{iSample,iModel} = bestParamTemp(idx,:);
            bestParamAllFm{iSample,iModel} = bestParamTemp;
            negLogLikeAllFm{iSample,iModel} = fvalTemp;
%             % Finishing timer and printing elapsed time
%             fprintf('fmincon elapsed time (days hours:minutes:seconds) %s \n\n',...
%                 datestr(etime(clock,cTemp)/86400,'dd HH:MM:SS'));
%             fprintf('Minimum negative log-likelihood: %.2f \n',fvalTemp(idx));
            
            %% Evaluating model fit
            
            % Save all simulation data and plot
            output = 'full';
            [negLogLike{iSample,iModel},bic{iSample,iModel},mdlEval{iSample,iModel}] = ...
                bci.fitModelAcrossRel(bestParamFm{iSample,iModel},paramNames,bciData,responseLoc,actModel,decisionFun,output);
            
            % Compute predicted distributions
            predDistr{iSample,iModel} = cat(1,mdlEval{iSample,iModel}.conditions);
            predDistr{iSample,iModel}.relV = cat(1,mdlEval{iSample,iModel}.relV);
            [predDistr{iSample,iModel}.f_predV,predDistr{iSample,iModel}.x_predV] = ...
                cellfun(funKs,{mdlEval{iSample,iModel}.sV_resp}','UniformOutput',false);
            [predDistr{iSample,iModel}.f_predA,predDistr{iSample,iModel}.x_predA] = ...
                cellfun(funKs,{mdlEval{iSample,iModel}.sA_resp}','UniformOutput',false);
            predDistr{iSample,iModel} = sortrows(predDistr{iSample,iModel},{'locV','locA','relV'});
            
            % Removing unnecessary fields from bciSimulations
            fieldsToRemove = {'sA_resp','sV_resp','sV_hat','sA_hat','s_hat_common',...
                'sV_hat_indep','sA_hat_indep'};
            for i = 1:numel(fieldsToRemove)
                if isfield(mdlEval{iSample,iModel},fieldsToRemove{i});
                    mdlEval{iSample,iModel} = rmfield(mdlEval{iSample,iModel},fieldsToRemove{i});
                end
            end
            
            % First fit null model
            negLogLike_null = bci.fitModelAcrossRel(bestParamFm{iSample,iModel},...
                                                    paramNames,bciData,responseLoc,'null',...
                                                    decisionFun,'logLike');
            % Compute R2 as well for the discrete case based on Nagelkerke (1991)
            n = size(bciData,1);
            R2{iSample,iModel} = 1-exp((-2/n)*(negLogLike_null-negLogLike{iSample,iModel}));
            maxR2 = 1-exp((2/n)*(-negLogLike_null));
            R2{iSample,iModel} = R2{iSample,iModel}/maxR2;
            
        end % end of loop over models
        
        respDistr{iSample} = varfun(wrapKs,bciData,...
            'InputVariables',{'respA','respV'},'GroupingVariables',{'locV','locA','relV'});
        respDistr{iSample}.f_respA = respDistr{iSample}.Fun_respA(:,1);
        respDistr{iSample}.x_respA = respDistr{iSample}.Fun_respA(:,2);
        respDistr{iSample}.f_respV = respDistr{iSample}.Fun_respV(:,1);
        respDistr{iSample}.x_respV = respDistr{iSample}.Fun_respV(:,2);
        respDistr{iSample}.Fun_respA = [];
        respDistr{iSample}.Fun_respV = [];
        respDistr{iSample}.GroupCount = [];
        respDistr{iSample}.Properties.RowNames = {};
        
        % Advancing progress monitor
        if progrMonitor, parfor_progress; end
        
    end % end of loop over samples
    
    % Closing progress monitor
    if progrMonitor, parfor_progress(0); end    
    
    % Compute predicted label variance
    plVar = num2cell(var(predLabel)');
        
    % Saving subject specific bci simulations
    dataFileMatchString = fileMatchStr;
    savePath = fullfile(workDir,['bci_simul_MVPA_',subID,'_',datestr(now,'yymmddHHMMSS'),'.mat']);
    save(savePath,'mdlEval','bestParamAllFm','bestParamFm',...
         'bestParamAllGs','bic','dataFileMatchString','negLogLike',...
         'negLogLikeAllGs','negLogLikeAllFm','models',...
         'neurDataSelection','paramCombinationsGs','parameterNames',...
         'plVar','predDistr','respDistr','R2','timePoints','-v7.3');
    if exist('isOutl','var')
        save(savePath,'isOutl,','nOutl','-append');
    end
    % Printing elapsed time
    cTemp = clock;
    fprintf('Elapsed time (days hours:minutes:seconds) %s \n\n',...
        datestr(etime(cTemp,cStart)/86400,'dd HH:MM:SS'));
end % end of loop over subjects

end
