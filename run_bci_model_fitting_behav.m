%% Preparations
clc;
clear all;

% Checking setup
[~,setupID] = system('hostname');
setupID = regexp(setupID,'[a-zA-Z_0-9-]*','match');
setupID = setupID{1};

if ~isempty(regexp(setupID,'^bb.+','once'))
    isServer = true;
    progrMonitor = false;
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

plotGridSerchFig = false;

%% Loading data
% subjList = {'109'};
subjList = {'108','109','110','111','112','113','116','118','119','120','121','122','123'};
expStage = 'final';

for iSubj = 1:numel(subjList)
    subID = subjList{iSubj};
    behavFilePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),...
        ['bci_data_BEHAV_',subID,'.mat']);
    if exist(behavFilePath,'file')
        bciData = load(behavFilePath);
        bciData = bciData.bciData;
    else
        warning('Missing behavioural data, returning');
        return;
    end
    
    %% Perform gridsearch on a grid of parameters
    % Response locations
    responseLoc = unique(bciData.locV);
    
    my_actions = 'logLike';  % only loglike compute
    
    % Preparing the parameter space
    % grid resolution
    nStep = 5;
    p_common = linspace(0.1,0.7,nStep);  %list of places to search for first parameter
    sigV1 = linspace(0.1,10,nStep);
    sigV2 = linspace(0.1,10,nStep);
    sigA = linspace(1,15,nStep);
    sigP = linspace(1,50,nStep);
    % These parameters can be used as well, but not included in this analysis
    % kernelWidth = 1; 'kW'
    % muP = 0; % mean of spatial prior
    
    % Parameter names
    parameterNamesBCI = {'p_common','sigP','sigA','sigV1','sigV2'};
    gridVectorsBCI = {p_common,sigP,sigA,sigV1,sigV2};
    nParameters = numel(gridVectorsBCI);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectorsBCI{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinationsBCI = cat(2,coords{:});
    
    % Generating parameter combinations for fusion model
    parameterNamesFus = {'sigP','sigA','sigV1','sigV2'};
    gridVectorsFus = {sigP,sigA,sigV1,sigV2};
    nParameters = numel(gridVectorsFus);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectorsFus{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinationsFus = cat(2,coords{:});
    
    % Task relevance model
    parameterNamesTaskRel = {'sigP','sigA','sigV1','sigV2'};
    paramCombinationsTaskRel = paramCombinationsFus;
    
    % Segregation visual model
    parameterNamesSegV = {'sigP','sigV1','sigV2'};
    gridVectorsSegV = {sigP,sigV1,sigV2};
    nParameters = numel(gridVectorsSegV);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectorsSegV{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinationsSegV = cat(2,coords{:});
    
    % Segregation visual model
    parameterNamesSegA = {'sigP','sigA'};
    gridVectorsSegA = {sigP,sigA};
    nParameters = numel(gridVectorsSegA);
    % Full factorial expansion of the specified parameter vectors
    coords = cell(1,nParameters);
    [coords{:}] = ndgrid(gridVectorsSegA{:});
    coords = cellfun(@(x) x(:),coords,'UniformOutput',false);
    % Matrix of all possible parameter combinations: nRows = number of
    % combinations, columns correspond to parameters
    paramCombinationsSegA = cat(2,coords{:});
    
    % Models to be fitted
    models = {'bci','bciSel','bciMatch','fus','taskRel'};
    % models = {'bci'};
    nModel = numel(models);
    
    % Pre-allocating arrays for data collection
    [negLogLike,negLogLikeGs,negLogLikeFm,bic,...
     mdlEval,bestParamGs,bestParamFm,R2] = deal(cell(1,nModel));
    paramCombinationsGs = {paramCombinationsBCI,paramCombinationsBCI,...
                        paramCombinationsBCI,paramCombinationsFus,...
                        paramCombinationsSegA,paramCombinationsSegV,...
                        paramCombinationsTaskRel};
    parameterNamesAll = {parameterNamesBCI,parameterNamesBCI,...
                        parameterNamesBCI,parameterNamesFus,...
                        parameterNamesSegA,parameterNamesSegV,...
                        parameterNamesTaskRel};
    
     for iModel = 1:nModel
        
        % Choosing model to fit
        actModel = models{iModel};
        
        switch actModel
          case 'bci'
            paramCombinations = paramCombinationsBCI;
            paramNames = parameterNamesBCI;
            % Decision function model averaging (1)
            % model selection (2)
            % probability matching (3)
            decisionFun = 1;
          case 'bciSel'
            paramCombinations = paramCombinationsBCI;
            paramNames = parameterNamesBCI;
            decisionFun = 2;
          case 'bciMatch'
            paramCombinations = paramCombinationsBCI;
            paramNames = parameterNamesBCI;
            decisionFun = 3;
          case 'fus'
            paramCombinations = paramCombinationsFus;
            paramNames = parameterNamesFus;
            decisionFun = 1;
          case 'segA'
            paramCombinations = paramCombinationsSegA;
            paramNames = parameterNamesSegA;
            decisionFun = 1;
          case 'segV'
            paramCombinations = paramCombinationsSegV;
            paramNames = parameterNamesSegV;
            decisionFun = 1;
          case 'taskRel'
            paramCombinations = paramCombinationsTaskRel;
            paramNames = parameterNamesTaskRel;
            decisionFun = 1;
        end
        
        % Pre-allocating array for data collection
        negLogLike_temp = NaN(size(paramCombinations,1),1);
        
        % Starting the timer and printing details
        cStart = clock;
        ds = datestr(now);
        printFnTitle(80,'run_bci_model_fitting',ds)
        fprintf('Performing gridsearch... \n');
        if progrMonitor, parfor_progress(size(paramCombinations,1)); end
        % Performing grid search on the specified parameter space
        parfor i = 1:size(paramCombinations,1)
            [negLogLike_temp(i)] = bci.fitModelAcrossRel(paramCombinations(i,:),...
                                                     paramNames,bciData,...
                                                     responseLoc,actModel,decisionFun,my_actions);
            if progrMonitor, parfor_progress; end
        end
        if progrMonitor, parfor_progress(0); end
        
        % Finding the minimum of the log-likelihood values and the corresponding
        % parameter combination
        [negLogLike_min,idx] = min(negLogLike_temp);
        startParFm = paramCombinations(idx,:);
        bestParamGs{iModel} = startParFm;
        
        if plotGridSerchFig
            figure;
            plot(negLogLike_temp,'k'); hold on; 
            plot(idx,negLogLike_temp(idx),'r*');
        end
        % Saving all log-likelihoods and parameter combinations
        negLogLikeGs{iModel} = negLogLike_temp;
        
        % Printing elapsed time
        cTemp = clock;
        fprintf('Gridsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
                datestr(etime(cTemp,cStart)/86400,'dd HH:MM:SS'));
        fprintf('Minimum negative log-likelihood: %.2f \n',negLogLike_min);
        
        %% Perform fminsearch to refine the best parameters
        % Take the best parameters found in the gridsearch and use them as starting
        % point in the fminsearch
        fprintf('Performing fminsearch... \n');
        % Options for fminsearch
        opts = optimset('fminsearch');
        opts.Display = 'final';
        opts.MaxFunEvals = 4000;
        opts.TolFun = 1.e-12;
        opts.MaxIter = 1000;
        
        % set upper and lower bounds on parameters
        switch actModel
          case 'bci'
            % 'p_common','sigP','sigA','sigV1','sigV2'
            LB = [0  0.001   0.001   0.001   0.001];
            UB = [1  inf     inf     inf     inf  ];
          case 'bciSel'
            % 'p_common','sigP','sigA','sigV1','sigV2'
            LB = [0  0.001   0.001   0.001   0.001];
            UB = [1  inf     inf     inf     inf  ];
          case 'bciMatch'
            % 'p_common','sigP','sigA','sigV1','sigV2'
            LB = [0  0.001   0.001   0.001   0.001];
            UB = [1  inf     inf     inf     inf  ];
          case 'fus'
            % 'sigP','sigA','sigV1','sigV2'
            LB = [0.001  0.001   0.001   0.001];
            UB = [inf    inf     inf     inf] ;
          case 'taskRel'
            % 'sigP','sigA','sigV1','sigV2'
            LB = [0.001  0.001   0.001   0.001];
            UB = [inf    inf     inf     inf];
          case 'segA'
            % 'sigP','sigA'
            LB = [0.001  0.001];
            UB = [inf    inf] ;
          case 'segV'
            % 'sigP','sigV1','sigV2'
            LB = [0.001  0.001   0.001];
            UB = [inf    inf     inf ] ;
        end
        
        % Creating anonymous function for input to fminsearchbnd
        fun = @(param) bci.fitModelAcrossRel(param,paramNames,bciData,...
                                             responseLoc,actModel,decisionFun,my_actions);
        [bestParamFm{iModel},negLogLikeFm{iModel}] = fminsearchbnd(fun,startParFm,LB,UB,opts);
        
        % Finishing timer and printing elapsed time
        fprintf('fminsearch elapsed time (days hours:minutes:seconds) %s \n\n',...
                datestr(etime(clock,cTemp)/86400,'dd HH:MM:SS'));
        
        %% Evaluating model fit
        
        % Save all simulation data and plot
        my_actions = 'full';
        
        % [logLike,bciSimulations] = bci.fitModelAcrossRel(bestParamFm,paramNames,...
                                                         % bciData,responseLoc,actModel,decisionFun,my_actions);
        [negLogLike{iModel},bic{iModel},mdlEval{iModel}] = ...
                bci.fitModelAcrossRel(bestParamFm{iModel},paramNames,bciData,...
                                      responseLoc,actModel,decisionFun,my_actions);
        
        % First fit null model
        negLogLike_null = bci.fitModelAcrossRel(bestParamFm{1},parameterNamesAll{1},bciData,...
                                                responseLoc,'null',decisionFun,'logLike');
        % Compute R2 as well for the discrete case based on Nagelkerke (1991)
        n = size(bciData,1);
        R2{iModel} = 1-exp((-2/n)*(negLogLike_null-negLogLike{iModel}));
        maxR2 = 1-exp((2/n)*(-negLogLike_null));
        R2{iModel} = R2{iModel}/maxR2;
     end % end of loop over models
        
    % Saving subject specific bci simulations
    fprintf('\n\nSaving data...\n\n');
    savePath = fullfile(DEC_2_setupdir(expStage,'anal_behav_sub',subID),...
        ['bci_simul_BEHAV_',subID,'_',datestr(now,'yymmddHHMMSS'),'.mat']);
    save(savePath,'mdlEval','bestParamFm',...
        'bestParamGs','bic','negLogLike',...
        'negLogLikeGs','negLogLikeFm','models',...
        'paramCombinationsGs','parameterNamesAll','R2','-v7.3');
    
end % End of loop over subjects


