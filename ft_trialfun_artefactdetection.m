function [trl,event] = ft_trialfun_artefactdetection(cfg)

% Finding the actual trials in the data based on the specified event values
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
stimOnsetEvSamples = evSamples(ismember(evValues,trStim));

% Removing sessions which are to be excluded
sessionToExclude = cfg.trialdef.fileSpec.exclude{:};
if ~isnan(sessionToExclude)
    for i = 1:numel(sessionToExclude)
        % Find and remove trials which belong to the session(s) to be
        % removed. 
        if sessionToExclude(i) ~= nSessions
            startExcluded = seStartEvSamples(sessionToExclude(i));
            endExcluded = seStartEvSamples(sessionToExclude(i)+1);
            stimOnsetEvSamples(stimOnsetEvSamples >= startExcluded & stimOnsetEvSamples < endExcluded) = [];
        else
            startExcluded = seStartEvSamples(sessionToExclude(i));
            stimOnsetEvSamples(stimOnsetEvSamples >= startExcluded) = [];
        end
    end
    % Removing sessions to be excluded and updating number of sessions
    seStartEvSamples(sessionToExclude) = [];
    nSessions = numel(seStartEvSamples);
end

se = NaN(nSessions,2);

% Finding session start and end samples. 
if nSessions == 1
    
    se(1) = seStartEvSamples;
    se(2) = stimOnsetEvSamples(end)+round(cfg.trialdef.poststim*hdr.Fs);
    
else
    
    se(:,1) = seStartEvSamples;
    
    for i = 1:nSessions
        
        % Find the last stimulus onset in the given session
        if i == nSessions
            lastStimInSessionSample = stimOnsetEvSamples(end);
        else
            lastStimInSessionSample  = max(stimOnsetEvSamples(stimOnsetEvSamples < se(i+1,1)));
        end
        se(i,2) = lastStimInSessionSample+round(cfg.trialdef.poststim*hdr.Fs);
        
    end

end

% Generating fake trials of length ~cfg.trialdef.trllength just for the 
% efficient artefact reviewing.
trlLengthSampl = cfg.trialdef.faketrllength*hdr.Fs;

% Number of artificial trials for each session. 
nTrialArtfPerSession = ceil((se(:,2)-se(:,1))/trlLengthSampl);

nTrialArtf = sum(nTrialArtfPerSession);

trl = zeros(nTrialArtf,3);

iTrial = 1;
for iSession = 1:nSessions
    
    startIdx = se(iSession,1);
    
    for iTrialArtfPerSession = 1:nTrialArtfPerSession(iSession)
        
        % Marking the beginning and ending samples of trials for artefact
        % reviewing. The offset remains always zero.
        trl(iTrial,1) = startIdx;
        trl(iTrial,2) = startIdx+(trlLengthSampl-1);
        startIdx = startIdx+trlLengthSampl;
        iTrial = iTrial+1;
        
    end
end

end