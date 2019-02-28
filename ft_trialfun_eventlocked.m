function [trl,event] = ft_trialfun_eventlocked(cfg)

hdr = ft_read_header(cfg.headerfile);
event = ft_read_event(cfg.headerfile);
evValues = {event.value}';
evTypes = {event.type}';
evSamples = [event.sample]';

isReqType = strcmp(evTypes,cfg.trialdef.eventtype);
evValues = evValues(isReqType);
evSamples = evSamples(isReqType);

% trigger values
trStim = cfg.trialdef.trigdef.trig_eeg(cfg.trialdef.trigdef.type == 'stim');
trResp = cfg.trialdef.trigdef.trig_eeg(cfg.trialdef.trigdef.type == 'resp');
trSessionStart = cfg.trialdef.trigdef.trig_eeg(cfg.trialdef.trigdef.type == 'sestart');
trTrialStart = cfg.trialdef.trigdef.trig_eeg(cfg.trialdef.trigdef.type == 'trstart');

% Correcting for trigger delay with respect to visual onset
temp = evSamples(ismember(evValues,[trStim',trTrialStart',trSessionStart']));
temp = temp+(cfg.trialdef.trig_visonset_corr*hdr.Fs);
evSamples(ismember(evValues,[trStim',trTrialStart',trSessionStart'])) = temp;

% all event values must be strings
nonStrIdx = find(~cellfun(@isstr,evValues));
evValues(nonStrIdx) = repmat({''},numel(nonStrIdx),1);

% Checking the number of sessions
nSessions = sum(ismember(evValues,trSessionStart));
if nSessions == 0
    error('No session start triggers were found!');
elseif nSessions ~= cfg.trialdef.fileSpec.nSessionsInFile
    error('Number of sessions does not match the expected!');
end
seStartEvSamples = evSamples(ismember(evValues,trSessionStart));
seStartEvValues = evValues(ismember(evValues,trSessionStart));
targEvSamples = evSamples(ismember(evValues,trStim));
targEvValues = evValues(ismember(evValues,trStim));

% Removing sessions which are to be excluded
sessionToExclude = cfg.trialdef.fileSpec.exclude{:};
if ~isnan(sessionToExclude)
    for i = 1:numel(sessionToExclude)
        % Find and remove trials which belong to the session(s) to be
        % removed. 
        if sessionToExclude(i) ~= nSessions
            startExcluded = seStartEvSamples(sessionToExclude(i));
            endExcluded = seStartEvSamples(sessionToExclude(i)+1);
            targEvValues(targEvSamples >= startExcluded & targEvSamples < endExcluded) = [];
            targEvSamples(targEvSamples >= startExcluded & targEvSamples < endExcluded) = [];
        else
            startExcluded = seStartEvSamples(sessionToExclude(i));
            targEvValues(targEvSamples >= startExcluded) = [];
            targEvSamples(targEvSamples >= startExcluded) = [];
        end
    end
    % Removing sessions to be excluded and updating number of sessions
    seStartEvSamples(sessionToExclude) = [];
    nSessions = numel(seStartEvSamples);
end

nTrials = numel(targEvSamples);

[begSamples,endSamples] = deal(NaN(nTrials,1));
[session,iTrialInSession,cond,hand,resp,respTime] = deal(NaN(nTrials,1));

actSessionInFile = 1;
actTrialInSession = 1;
if seStartEvValues{actSessionInFile} == trSessionStart(1)
    actHand = 1;
elseif seStartEvValues{actSessionInFile} == trSessionStart(2)
    actHand = 2;
end

for i = 1:nTrials
    
    % Find actual target event sample
    actTargEvSampl = targEvSamples(i);
    % Updating actual session and actual trial in session
    if actSessionInFile ~= find(actTargEvSampl > seStartEvSamples,1,'last')
        actSessionInFile = find(actTargEvSampl > seStartEvSamples,1,'last');
        if seStartEvValues{actSessionInFile} == trSessionStart(1)
            actHand = 1;
        elseif seStartEvValues{actSessionInFile} == trSessionStart(2)
            actHand = 2;
        end
        actTrialInSession = 1;
    end
    % Find the index of the sample corresponding to the target event of
    % the trial.
    actTargEvIdx = find(evSamples == actTargEvSampl);
    % Find the index of the sample corresponding to the response to the
    % trial. 
    noResponse = false;
    actRespEvIdx = [];   
    j = 1; 
    while 1
        
        anEvVal = evValues(actTargEvIdx+j);
        if any(ismember(trResp,anEvVal))
            actRespEvIdx = actTargEvIdx+j;
            break;            
        elseif any(ismember([trTrialStart',trSessionStart'],anEvVal))
            noResponse = true;
            break;
        elseif actTargEvIdx+j >= size(evValues,1)
            noResponse = true;
            break;
        end
        
        j = j+1;
    end
    
    % Beginning and ending samples of trials
    begSamples(i) = actTargEvSampl-round(cfg.trialdef.prestim*hdr.Fs);
    endSamples(i) = actTargEvSampl+round(cfg.trialdef.poststim*hdr.Fs);
    % Session number
    session(i) = cfg.trialdef.fileSpec.iSessionOverall{:}(actSessionInFile);
    % Serial number of trial in the session
    iTrialInSession(i) = actTrialInSession;
    % Stimulus condition
    cond(i) = cfg.trialdef.trigdef.cond(cfg.trialdef.trigdef.trig_eeg == targEvValues{i});
    % Hand to respond
    hand(i) = actHand;
    % Response and respTime
    if noResponse
        resp(i) = NaN;
        respTime(i) = NaN;
    else
        resp(i) = cfg.trialdef.trigdef.loc(cfg.trialdef.trigdef.trig_eeg == evValues(actRespEvIdx));
        respTime(i) = (evSamples(actRespEvIdx)-actTargEvSampl)/hdr.Fs;
    end
    
    actTrialInSession = actTrialInSession + 1;
end

% Creating the offset array
offset = -round(cfg.trialdef.prestim*hdr.Fs)*ones(size(begSamples));

trl = [begSamples,endSamples,offset,session,iTrialInSession,cond,hand,resp,respTime];

end