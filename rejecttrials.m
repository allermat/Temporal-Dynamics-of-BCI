function badTrials = rejecttrials(eegData,behavData,eyeData,artfData,trigDef,winOfInt,minRespTime,saccadeThresholds,fixationCriteria,mode)
% Checks if trials meet various criteria
% 
% DETAILS:
%   Checks each trial found in the eegData whether: 
%       (1) the response is missing
%       (2) the response is early
%       (3) the response was made by wrong hand    
%       (4) there is an EEG artefact
%       (5) there is eyetracker data present for the trial
%       (6) there is an eyeblink
%       (7) there is a saccade
%       (8) the fixation is at the correct location
%       
% INPUT:
%   eegData: epoched fieldtrip data file
%   behavData: table of behavioural data (containing the whole dataset)
%   eyeData: structure array of eyetracker data (containing the whole 
%            dataset)
%   artfData: structure array of detected EEG artefacts (corresponding to
%       the fieldtrip data file)
%   trigDef: table of trigger code definitions
%   winOfInt: 1 x 2 vector with the start and end time of time window of 
%       interest with respect to the target event onset
%   minRespTime: scalar, minimum response time in seconds
%   saccadeThresholds: 1 x 3 vector of the threshold values for defining 
%       the saccades
%       thresholds(1) = minimum amplitude in degrees
%       thresholds(2) = minimum peak velocity in degrees/seconds
%       thresholds(3) = minimum duration in seconds
%   fixationCriteria: 1 x 3 vector of the criteria for accepting a fixation
%       event
%       criteria(1) = fixation cross x coordinate (in pixels)
%       criteria(2) = fixation cross y coordinate (in pixels)
%       criteria(3) = radius of the accepted area around the
%                     fixation cross (in pixels)
%   mode: 'normal', 'prestim'
% 
% OUTPUT: 
%   badTrials: vector of number of trials x 1 with the following codes
%       0 - good
%       1 - no response
%       2 - early response
%       3 - wrong hand
%       4 - EEG artefact
%       5 - missing eyetracker data
%       6 - eyeblink
%       7 - saccade
%       8 - wrong fixation location
% 

%% Parsing input
p = inputParser;

% Input checking functions
checkFtData = @(x) ft_datatype(x,'raw');
checkBehavData = @(x) istable(x) && ...
    all(ismember({'iSessionOverall','iTrialInSession','respTime'},x.Properties.VariableNames));
checkEyeData = @(x) isstruct(x) && all(ismember({'session','event','Fs'},fieldnames(x)));
checkArtfData = @(x) isstruct(x) && all(ismember(fieldnames(x),{'artefact_muscle','artefact_visual'}));
checkTrigDef = @(x) istable(x) && ...
    all(ismember({'trig_eeg','trig_eye','type','cond','loc'},x.Properties.VariableNames));
validModes = {'normal','prestim'};

% Defining input
addRequired(p,'eegData',checkFtData);
addRequired(p,'behavData',checkBehavData);
addRequired(p,'eyeData',checkEyeData);
addRequired(p,'artfData',checkArtfData);
addRequired(p,'trigDef',checkTrigDef);
addRequired(p,'winOfInt',@(x)validateattributes(x,{'numeric'},{'size',[1,2],'finite','nonnan'}));
addRequired(p,'minRespTime',@(x)validateattributes(x,{'numeric'},{'scalar','positive','finite','nonnan'}));
addRequired(p,'saccadeThresholds',@(x)validateattributes(x,{'numeric'},{'size',[1,3],'positive','finite','nonnan'}));
addRequired(p,'fixationCriteria',@(x)validateattributes(x,{'numeric'},{'size',[1,3],'positive','finite','nonnan'}));
addRequired(p,'mode',@(x)any(validatestring(x,validModes)));

% Parsing inputs
parse(p,eegData,behavData,eyeData,artfData,trigDef,winOfInt,minRespTime,...
    saccadeThresholds,fixationCriteria,mode);

% Assigning input to variables
eegData = p.Results.eegData;
behavData = p.Results.behavData;
eyeData = p.Results.eyeData;
artfData = p.Results.artfData;
trigDef = p.Results.trigDef;
winOfInt = p.Results.winOfInt;
minRespTime = p.Results.minRespTime;
saccadeThresholds = p.Results.saccadeThresholds;
fixationCriteria = p.Results.fixationCriteria;
mode = p.Results.mode;

%% Function body

% eeg data properties
nTrials = size(eegData.trial,2);
trialInfo = eegData.trialinfo;
actSession = trialInfo(1,1);

% Array for collecting rejection info
badTrials = zeros(nTrials,1);

for iTrialInFile = 1:nTrials
    
    % Updating the actual session and related variables if necessary
    if actSession ~= trialInfo(iTrialInFile,1)
        actSession = trialInfo(iTrialInFile,1);
    end
    
    actTrialInSession = trialInfo(iTrialInFile,2);
    
    %% Checking response
    if strcmp(mode,'normal')
        % Is there a response?
        if isnan(behavData.resp(behavData.iSessionOverall == actSession & behavData.iTrialInSession == actTrialInSession))
            badTrials(iTrialInFile) = 1;
            continue;
        end
        
        % Is the response too early?
        if behavData.respTime(behavData.iSessionOverall == actSession & behavData.iTrialInSession == actTrialInSession) < minRespTime
            badTrials(iTrialInFile) = 2;
            continue;
        end
        
        % Is the response too early?
        if behavData.wrongHand(behavData.iSessionOverall == actSession & behavData.iTrialInSession == actTrialInSession)
            badTrials(iTrialInFile) = 3;
            continue;
        end
    end
    %% Checking EEG artefacts
    % Is there an EEG artefact?
    if checkartefacts(eegData.sampleinfo(iTrialInFile,:),eegData.time{iTrialInFile},eegData.fsample,winOfInt,artfData)
        badTrials(iTrialInFile) = 4;
        continue;
    end
    
    %% Checking eye events
    % Finding eye events corresponding the current trial
    eyeEventsActTrial = selecteyeevents(eyeData(actSession),actTrialInSession,winOfInt,trigDef);
    if isempty(eyeEventsActTrial)
        warning('No eyetracking data for session %d, trial %d',actSession,actTrialInSession);
        badTrials(iTrialInFile) = 5;
        continue;
    end
    
    % Is there a blink? 
    if any(ismember({eyeEventsActTrial.type},'Blink'))
        badTrials(iTrialInFile) = 6;
        continue;
    end
    
    % Is there a saccade? 
    if any(ismember({eyeEventsActTrial.type},'Saccade'))
        
%         % Using saccades as they are defined by the EyeLink software
%         badTrials(iTrialInFile) = 7;
%         continue;
        
        % Further filtering saccades
        saccades = eyeEventsActTrial(ismember({eyeEventsActTrial.type},'Saccade'));
        % Are any of the saccades above threshold?
        if checksaccades(saccades,saccadeThresholds,eyeData(actSession).Fs)
            badTrials(iTrialInFile) = 7;
            continue;
        end
        
    end
    
    % Is fixation average location correct?
    if all(~ismember({eyeEventsActTrial.type},'Fixation'))
        badTrials(iTrialInFile) = 8;
        continue;
    else
        fixations = eyeEventsActTrial(ismember({eyeEventsActTrial.type},'Fixation'));
        % Are any of the saccades above threshold?
        if checkfixations(fixations,fixationCriteria)
            badTrials(iTrialInFile) = 8;
        end
        
    end
    
end

end


function foundArtefact = checkartefacts(sampleInfo,time,Fs,winOfInt,artfData)
% Checks whether there is any artefact within the trial
% 
% INPUT: 
%   sampleInfo: 1 x 2 vector of the start and end samples of the trial of
%       interest
%   time: 1 x N vector of time values (in seconds) corresponding to the 
%       samples of the trial of interest, where N is the number of samples.
%   Fs: sampling frequency
%   winOfInt: 1 x 2 vector with the start and end time of time window of 
%       interest with respect to the target event onset
%   artfData: structure array of artefact definitions. Each field of the
%       array is an N x 2 matrix of artefact start and end samples, where N
%       is the number of artefacts. Different fields contain different
%       artefact types. 
% 
% OUTPUT:
%   foundArtefact: true if there is an artefact within the window of
%       interest of the trial
% 

foundArtefact = false;

targEventOnsetSample = sampleInfo(1)-(time(1)*Fs);
startSampleWinOfInt = targEventOnsetSample+(winOfInt(1)*Fs);
endSampleWinOfInt = targEventOnsetSample+(winOfInt(2)*Fs);

artfTypes = fieldnames(artfData);

for i = 1:size(artfTypes,1)
    
    actArtfData = artfData.(artfTypes{i});
    
    % Skip if there are no artefacts
    if isempty(actArtfData)
        break;
    end
    
    % Find artefacts which:
    % begin within the window of interest
    crit1 = actArtfData(:,1) >= startSampleWinOfInt & actArtfData(:,1) <= endSampleWinOfInt;
    % end within the winow of interest
    crit2 = actArtfData(:,2) >= startSampleWinOfInt & actArtfData(:,2) <= endSampleWinOfInt;
    % cover the whole window of interest
    crit3 = actArtfData(:,1) <= startSampleWinOfInt & actArtfData(:,2) >= endSampleWinOfInt;
    
    % If there is an artefact which fall under either criteria, mark it and
    % break. 
    foundArtefact = any(crit1 | crit2 | crit3);
    if foundArtefact
        break;
    end
    
end

end


function selectedEvents = selecteyeevents(eyeData,actTrial,winOfInt,trigDef)
% Selects eyetracker events corresponding to a particular trial
% 
% INPUT: 
%   eyeData: structure array of eyetracker data of the session
%       corresponding to the trial of interest
%   actTrial: serial number of the trial of interest within its session
%   winOfInt: 1 x 2 vector with the start and end time of time window of 
%       interest with respect to the target event onset
%   trigDef: trigger definition table
% 
% OUTPUT:
%   selectedEvents: structure array of events corresponding to the trial of
%       interest. Empty structure if no events were found
% 

evValues = {eyeData.event.value}';
evTypes = {eyeData.event.type}';
evStartSamples = [eyeData.event.sample]';
evEndSamples = evStartSamples+[eyeData.event.duration]';
targEvStartSamples = evStartSamples(ismember(evTypes,'Stimulus') & ...
    ismember(evValues,trigDef.trig_eye(trigDef.type == 'stim')));

% If the trial is not found return an empty structure
if ~ismember(1:size(targEvStartSamples,1),actTrial)
    selectedEvents = struct([]);
    return;
end

% Find window of interest around target event
startSampleWinOfInt = targEvStartSamples(actTrial)+(winOfInt(1)*eyeData.Fs);
endSampleWinOfInt = targEvStartSamples(actTrial)+(winOfInt(2)*eyeData.Fs);

% Find events which:
% begin within the window of interest
crit1 = evStartSamples >= startSampleWinOfInt & evStartSamples <= endSampleWinOfInt;
% end within the winow of interest
crit2 = evEndSamples >= startSampleWinOfInt & evEndSamples <= endSampleWinOfInt;
% cover the whole window of interest
crit3 = evStartSamples <= startSampleWinOfInt & evEndSamples >= endSampleWinOfInt;
    
% Choose those events which fall under any of the above criteria
selectedEvents = eyeData.event(crit1 | crit2 | crit3);

end


function aboveThreshold = checksaccades(saccades,thresholds,Fs)
% Checks whether saccade parameters are above the given threshold
% 
% INPUT: 
%   saccades: structure array containing saccade events
%   thresholds: 1 x 3 vector of the threshold values for defining a saccade
%               thresholds(1) = minimum amplitude in degrees
%               thresholds(2) = minimum peak velocity in degrees/seconds
%               thresholds(3) = minimum duration in seconds
%   Fs: sampling frequency
% 
% OUTPUT:
%   aboveThreshold: true if any of the saccade events are above all three 
%       thresholds. 
% 

if any(~ismember({saccades.type},'Saccade'))
    error('Event of wrong type passed to checksaccades!');
end

temp = regexp({saccades.value},'([0-9.e+]+)\s+([0-9.e+]+)','tokens');
if isempty(temp)
    error('Wrong saccade value format!');
end
temp = [temp{:}]';
temp = vertcat(temp{:});
temp = cell2mat(cellfun(@str2num,temp,'UniformOutput',false));
ampl = temp(:,1);
pv = temp(:,2);
% Converting duration from samples to seconds
dur = [saccades.duration]'./Fs;

if ampl >= thresholds(1) % & pv >= thresholds(2) & dur >= thresholds(3))
    aboveThreshold = true;
else
    aboveThreshold = false;
end

end


function failedCriteria = checkfixations(fixations,criteria)
% Checks whether fixation parameters meet the given criteria
% 
% INPUT: 
%   fixations: structure array containing fixation events
%   criteria: vector of the criteria for accepting a fixation event
%               criteria(1) = fixation cross x coordinate (in pixels)
%               criteria(2) = fixation cross y coordinate (in pixels)
%               criteria(3) = radius of the accepted area around the
%                             fixation cross (in pixels)
% 
% OUTPUT:
%   failedCriteria: true if any of the fixation events fail to meet the
%       criteria
% 

if any(~ismember({fixations.type},'Fixation'))
    error('Event of wrong type passed to checkfixations!');
end

temp = regexp({fixations.value},'([0-9.e+]+)\s+([0-9.e+]+)','tokens');
if isempty(temp)
    error('Wrong saccade value format!');
end
temp = [temp{:}]';
temp = vertcat(temp{:});
temp = cell2mat(cellfun(@str2num,temp,'UniformOutput',false));
% This should produce an N x 2 matrix, where N is the number of fixation
% events and the two columns are the average x and y positions 
% respectively. 

if any(~isPointInCircle(temp,repmat(criteria,size(temp,1),1)))
    failedCriteria = true;
else
    failedCriteria = false;
end

end