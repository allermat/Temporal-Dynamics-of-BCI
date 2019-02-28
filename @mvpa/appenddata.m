function appenddata(I)
% Method for appending data to an existing MVPA dataset
%
% INPUT: 
%   I (struct): input settings. 
%       Required fields:
%           pathFileToReceive (string): path to the dataset which receives 
%               the data to be appended.
%           to_append (string): name of the data to be appended. Possible
%               values: 'prestim_alpha_data'
%       Optional fields:
%           pathFileToAppend (string): path to the dataset which contains 
%               the data to be appended.
% OUTPUT: -
%   The dataset will be saved to disk with the appended data

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;

validToAppend = {'bci_model_data','prestim_alpha_data','prestim_beta_data',...
    'prestim_theta_data','locAxorVlabel','finger','LvsRa'};

addParameter(p,'pathFileToReceive','',@(x) exist(x,'file'));
addParameter(p,'pathFileToAppend','',@(x) exist(x,'file'));
addParameter(p,'to_append','',@(x)any(validatestring(x,validToAppend)));

parse(p,I);

pathFileToReceive = p.Results.pathFileToReceive;
pathFileToAppend = p.Results.pathFileToAppend;
toAppend = p.Results.to_append;

if ismember(toAppend,{'bci_model_data','prestim_alpha_data',...
                      'prestim_beta_data','prestim_theta_data'})
    requiredVars = {'pathFileToReceive','pathFileToAppend','to_append'};
else
    requiredVars = {'pathFileToReceive','to_append'};
end

if any(ismember(requiredVars,p.UsingDefaults))
    error('All required parameters must be specified!');
end

% Open file to append data to using matfile  
dataFileToReceive = matfile(pathFileToReceive,'Writable',true);
dataFileInfoToReceive = dataFileToReceive.info;

% Preparing data to be appended
switch toAppend
    case 'bci_model_data'
        dataToReceiveVarNames = {'sA_hat_and_sV_hat','s_hat_comm','sV_hat_indep','sA_hat_indep'};
        if any(ismember(dataToReceiveVarNames,dataFileInfoToReceive.Properties.VariableNames))
            warning('The data to be appended is already part of the dataset, overwriting. ');
        end
        % Pre-allocating new variables
        for i = 1:numel(dataToReceiveVarNames)
            dataFileInfoToReceive.(dataToReceiveVarNames{i}) = NaN(size(dataFileInfoToReceive,1),1);
        end
        % Opening data file to append
        dataToAppend = load(pathFileToAppend);
        dataToAppend = dataToAppend.bciSimulations;
        % Extracting conditions
        conditions = cat(1,dataToAppend.conditions);
        conditions.relV = [dataToAppend.relV]';
        % Looping thorugh conditions and saving variables
        for i = 1:size(conditions,1)
            actCondIdx = ismember(dataFileInfoToReceive(:,{'locV','locA','relV'}),conditions(i,:),'rows');
            for j = 1:2
                % loop thorug task levels
                actCondIdxTask = actCondIdx & dataFileInfoToReceive.task == j;
                switch j
                    case 1, sHatStr = 'sA_hat';
                    case 2, sHatStr = 'sV_hat';
                end
                dataFileInfoToReceive.sA_hat_and_sV_hat(actCondIdxTask) = dataToAppend(i).(sHatStr);
                dataFileInfoToReceive.s_hat_comm(actCondIdxTask) = dataToAppend(i).s_hat_common;
                dataFileInfoToReceive.sV_hat_indep(actCondIdxTask) = dataToAppend(i).sV_hat_indep;
                dataFileInfoToReceive.sA_hat_indep(actCondIdxTask) = dataToAppend(i).sA_hat_indep;
            end
        end
        % Append new variables to the dataset
        dataFileToReceive.info = dataFileInfoToReceive;
    case 'finger'
        dataToReceiveVarNames = {'finger'};
        if any(ismember(dataToReceiveVarNames,dataFileInfoToReceive.Properties.VariableNames))
            warning('The data to be appended is already part of the dataset, overwriting. ');
        end
        % Pre-allocating new variables
        dataFileInfoToReceive.finger = ...
            NaN(size(dataFileInfoToReceive,1),1);
        u = unique(dataFileInfoToReceive(:,{'hand','resp'}),'rows');
        for i = 1:size(u,1)
            idx = ismember(dataFileInfoToReceive(:,{'hand','resp'}),u(i,:));
            dataFileInfoToReceive.finger(idx) = i;
        end
        if any(isnan(dataFileInfoToReceive.finger))
            warning('mpva:appenddata:missingValues', ...
                    'There are missing values in the ''finger'' variable');
        end
        % Append new variables to the dataset
        dataFileToReceive.info = dataFileInfoToReceive;
    case 'locAxorVlabel'
        dataToReceiveVarNames = {'locAxorV'};
        if any(ismember(dataToReceiveVarNames,dataFileInfoToReceive.Properties.VariableNames))
            warning('The data to be appended is already part of the dataset, overwriting. ');
        end
        % Pre-allocating new variables
        dataFileInfoToReceive.locAxorV = NaN(size(dataFileInfoToReceive,1),1);
        dataFileInfoToReceive.locAxorV(isnan(dataFileInfoToReceive.locA)) = ...
            dataFileInfoToReceive.locV(isnan(dataFileInfoToReceive.locA));
        dataFileInfoToReceive.locAxorV(isnan(dataFileInfoToReceive.locV)) = ...
            dataFileInfoToReceive.locA(isnan(dataFileInfoToReceive.locV));
        % Append new variables to the dataset
        dataFileToReceive.info = dataFileInfoToReceive;
    case 'LvsR'
        dataToReceiveVarNames = {'LvsRa','LvsRv'};
        if any(ismember(dataToReceiveVarNames,dataFileInfoToReceive.Properties.VariableNames))
            warning('The data to be appended is already part of the dataset, overwriting. ');
        end
        % Pre-allocating new variables
        dataFileInfoToReceive.LvsRa = NaN(size(dataFileInfoToReceive,1),1);
        dataFileInfoToReceive.LvsRa(dataFileInfoToReceive.locA < 0 ) = -1;
        dataFileInfoToReceive.LvsRa(dataFileInfoToReceive.locA > 0 ) = 1;
        dataFileInfoToReceive.LvsRv = NaN(size(dataFileInfoToReceive,1),1);
        dataFileInfoToReceive.LvsRv(dataFileInfoToReceive.locV < 0 ) = -1;
        dataFileInfoToReceive.LvsRv(dataFileInfoToReceive.locV > 0 ) = 1;
        % Append new variables to the dataset
        dataFileToReceive.info = dataFileInfoToReceive;
        
    case 'prestim_alpha_data'
        dataToReceiveVarName = 'psAlphaPowOcc';
        if ismember(dataToReceiveVarName,dataFileInfoToReceive.Properties.VariableNames)
            warning('The data to be appended is already part of the dataset, overwriting. ');
        end
        % Opening data file to append
        dataToAppend = load(pathFileToAppend);
        dataToAppend = dataToAppend.ftDataPow;
        % Extracting prestim alpha averaged over occipital electrodes
        isAlpha = dataToAppend.freq >= 7.5 & dataToAppend.freq <= 14.5;
        isOccElec = ismember(dataToAppend.label,{'PO3','POz','PO4','O1','Oz','O2'});
        temp = dataToAppend.powspctrm;
        temp = mean(temp(:,:,isAlpha),3);
        psAlphaOccPowByTrial = mean(temp(:,isOccElec),2);
        % Finding the trials in the receiving file which are also present
        % in the to be appened file
        tempTable = table(dataToAppend.trialinfo(:,1),dataToAppend.trialinfo(:,2),'VariableNames',{'session','iTrialSession'});
        [Lia,Locb] = ismember(dataFileInfoToReceive(:,{'session','iTrialSession'}),tempTable);
        % Appending the data
        dataFileInfoToReceive.(dataToReceiveVarName) = NaN(size(dataFileInfoToReceive,1),1);
        dataFileInfoToReceive.(dataToReceiveVarName)(Lia) = psAlphaOccPowByTrial(Locb(Locb ~= 0));
        dataFileToReceive.info = dataFileInfoToReceive;
        
    case 'prestim_beta_data'
        dataToReceiveVarName = 'psBetaPowTemporoPar';
        if ismember(dataToReceiveVarName,dataFileInfoToReceive.Properties.VariableNames)
            warning('The data to be appended is already part of the dataset, overwriting. ');
        end
        % Opening data file to append
        dataToAppend = load(pathFileToAppend);
        dataToAppend = dataToAppend.ftDataPow;
        % Extracting prestim beta averaged over occipital electrodes
        isBeta = dataToAppend.freq >= 15 & dataToAppend.freq <= 30;
        elecSelection = ismember(dataToAppend.label,...
            {'T7','C5','C3','C1','Cz','C2','C4','C6','T8',...
            'TP9','TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6','TP8','TP10',...
            'P7','P5','P3','P1','Pz','P2','P4','P6','P8'});
        temp = dataToAppend.powspctrm;
        temp = mean(temp(:,:,isBeta),3);
        psBetaPowByTrial = mean(temp(:,elecSelection),2);
        % Finding the trials in the receiving file which are also present
        % in the to be appened file
        tempTable = table(dataToAppend.trialinfo(:,1),dataToAppend.trialinfo(:,2),'VariableNames',{'session','iTrialSession'});
        [Lia,Locb] = ismember(dataFileInfoToReceive(:,{'session','iTrialSession'}),tempTable);
        % Appending the data
        dataFileInfoToReceive.(dataToReceiveVarName) = NaN(size(dataFileInfoToReceive,1),1);
        dataFileInfoToReceive.(dataToReceiveVarName)(Lia) = psBetaPowByTrial(Locb(Locb ~= 0));
        dataFileToReceive.info = dataFileInfoToReceive;
        
    case 'prestim_theta_data'
        dataToReceiveVarName = 'psThetaPowFrontal';
        if ismember(dataToReceiveVarName,dataFileInfoToReceive.Properties.VariableNames)
            warning('The data to be appended is already part of the dataset, overwriting. ');
        end
        % Opening data file to append
        dataToAppend = load(pathFileToAppend);
        dataToAppend = dataToAppend.ftDataPow;
        % Extracting prestim beta averaged over occipital electrodes
        isTheta = dataToAppend.freq >= 4 & dataToAppend.freq <= 7;
        elecSelection = ismember(dataToAppend.label,...
            {'AF','AF3','AF4','AF8',...
            'F7','F5','F3','F1','Fz','F2','F4','F6','F8',...
            'FT7','FC5','FC3','FC1','FC2','FC4','FC6','FT8'});
        temp = dataToAppend.powspctrm;
        temp = mean(temp(:,:,isTheta),3);
        psThetaPowByTrial = mean(temp(:,elecSelection),2);
        % Finding the trials in the receiving file which are also present
        % in the to be appened file
        tempTable = table(dataToAppend.trialinfo(:,1),dataToAppend.trialinfo(:,2),'VariableNames',{'session','iTrialSession'});
        [Lia,Locb] = ismember(dataFileInfoToReceive(:,{'session','iTrialSession'}),tempTable);
        % Appending the data
        dataFileInfoToReceive.(dataToReceiveVarName) = NaN(size(dataFileInfoToReceive,1),1);
        dataFileInfoToReceive.(dataToReceiveVarName)(Lia) = psThetaPowByTrial(Locb(Locb ~= 0));
        dataFileToReceive.info = dataFileInfoToReceive;
        
    otherwise
        error('This data type is not yet supported! ');
end

end
