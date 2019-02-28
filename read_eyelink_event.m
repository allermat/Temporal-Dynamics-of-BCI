function [event,hdr] = read_eyelink_event(fileName)
% Read MSG, ESACC, EBLINK, EFIX and INPUT events from eyelink .asc files
% 
% USAGE: 
%   [event,hdr] = read_eyelink_event(fileName)
% 
% DETAILS: 
%   This piece of code assumes that the recording contains a single block
%       of continuous recording without any break. If your recording 
%       contains multiple blocks of recordings with breaks intermingled you
%       might encounter errors and funny behaviour. 
%   ft_read_event does not read MSG, SACC and FIX events into the event
%       structure as the syntax of those events can vary massively. This
%       function reads these events but with only the below syntax (see
%       EyeLink manual for other syntaxes)
%   
%   MSG syntax: 
%       MSG <time> <message> 
%           (I assume, that <message> contains only digits)
%   
%   ESACC syntax: 
%       ESACC <eye> <stime> <etime> <dur> <sxp> <syp> <exp> <eyp> <ampl> <pv>
%   
%   EBLINK syntax: 
%       EBLINK <eye> <stime> <etime> <dur>
%   
%   EFIX syntax: 
%       EFIX <eye> <stime> <etime> <dur> <axp> <ayp> <aps>
%   
%   NOTATIONS
%       <eye>         which eye caused event ("L" or "R")
%       <time>        timestamp in milliseconds
%       <stime>       timestamp of first sample in milliseconds
%       <etime>       timestamp of last sample in milliseconds
%       <dur>         duration in milliseconds
%       <axp>, <ayp>  average X and Y position
%       <sxp>, <syp>  start X and Y position data
%       <exp>, <eyp>  end X and Y position data
%       <aps>         average pupil size (area or diameter)
%       <av>, <pv>    average, peak velocity (degrees/sec)
%       <ampl>        saccadic amplitude (degrees)
%       <xr>, <yr>    X and Y resolution (position units/degree)
% 
% INPUT: 
%   fileName: name of the file with extension (or full path)
% 
% OUTPUT:
%   event: structure array of events
%   hdr: header of the file as returned by ft_read_header
% 

%% Parsing input, additional checks
p = inputParser;

checkFileName = @(x) exist(x,'file');
addRequired(p,'fileName',checkFileName);

parse(p,fileName);

fileName = p.Results.fileName;


% Check for fieldtrip
if isempty(regexp(path,'fieldtrip','once'))
    error('fieldtrip toolbox was not found in path');
end
% Check extension;
[~,~,ext] = fileparts(fileName);
if ~strcmp(ext,'.asc')
    error('Wrong file type! ');
end

%% Reading raw data
hdr = ft_read_header(fileName);
% The built in function doesn't parse the lines of blink events, so we have
% to replace it with our own version with bilnks added. (See below)
hdr.orig = read_eyelink_asc_ext(fileName);
% Message events
msg = hdr.orig.msg;
% End of saccade events
sacc = hdr.orig.esacc;
% End of blink events
blink = hdr.orig.eblink;
% End of fixation events
fixx = hdr.orig.efix;
% Input events
inp = hdr.orig.input;
% Finding stimulus events amongst message events
msg = regexp(msg,'MSG\s\d+\s\d+','match','once');
msg = msg(~strcmp(msg,''));

% Initializing event array
nMsg = size(msg,1);
nSacc = size(sacc,1);
nBlink = size(blink,1);
nFix = size(fixx,1);
nInp = size(inp,1);
event = struct('type','','value','','sample',[],'duration',[],'offset',[],'timestamp',[]);
event(nMsg+nSacc+nBlink+nFix+nInp).type = '';

% Anonymous function to convert timestamp values to samples
tstamp2sample = @(x) round((x-hdr.FirstTimeStamp)./hdr.TimeStampPerSample) + 1;

% Extracting MSG details
% MSG <time> <message> 
if nMsg > 0
    temp = regexp(msg,'MSG\s+(\d+)\s+(\d+)','tokens');
    temp = [temp{:}]';
    temp = vertcat(temp{:});
    temp(:,1) = cellfun(@str2num,temp(:,1),'UniformOutput',false);
    % Converting timestamps to samples
    temp(:,3) = cellfun(tstamp2sample,temp(:,1),'UniformOutput',false);
    % Saving MSG event as Stimulus type
    startIdx = 1;
    endIdx = nMsg;
    [event(startIdx:endIdx).type] = deal('Stimulus');
    [event(startIdx:endIdx).sample] = temp{:,3};
    [event(startIdx:endIdx).value] = temp{:,2};
    [event(startIdx:endIdx).duration] = deal(1);
    [event(startIdx:endIdx).offset] = deal(0);
    [event(startIdx:endIdx).timestamp] = temp{:,1};
end

% Extracting ESACC details
% ESACC <eye> <stime> <etime> <dur> <sxp> <syp> <exp> <eyp> <ampl> <pv>
% For floating point values the exponential notation is also possible, so I
% included that in the regular expression. 
if nSacc > 0
    temp = regexp(sacc,...
        'ESACC\s+[LR]\s+(\d+)\s+(\d+)\s+[0-9.e+]+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+[0-9.e+-]+\s+([0-9.e+]+\s+[0-9.e+]+)','tokens');
    temp = [temp{:}]';
    temp = vertcat(temp{:});
    temp(:,1:2) = cellfun(@str2num,temp(:,1:2),'UniformOutput',false);
    % Converting timestamps to samples
    temp(:,4:5) = cellfun(tstamp2sample,temp(:,1:2),'UniformOutput',false);
    % Computing the duration in samples
    temp(:,6) = cellfun(@minus,temp(:,5),temp(:,4),'UniformOutput',false);
    % Saving saccade events
    % The amlitude and peak velocity are saved as the vent value
    startIdx = nMsg+1;
    endIdx = nMsg+nSacc;
    [event(startIdx:endIdx).type] = deal('Saccade');
    [event(startIdx:endIdx).sample] = temp{:,4};
    [event(startIdx:endIdx).value] = temp{:,3};
    [event(startIdx:endIdx).duration] = temp{:,6};
    [event(startIdx:endIdx).offset] = deal(0);
    [event(startIdx:endIdx).timestamp] = temp{:,1};
end

% Extracting EBLINK details
% EBLINK <eye> <stime> <etime> <dur>
if nBlink > 0
    temp = regexp(blink,'EBLINK\s+[LR]\s+(\d+)\s+(\d+)\s+[0-9.e+]','tokens');
    temp = [temp{:}]';
    temp = vertcat(temp{:});
    temp(:,1:2) = cellfun(@str2num,temp(:,1:2),'UniformOutput',false);
    % Converting timestamps to samples
    temp(:,3:4) = cellfun(tstamp2sample,temp(:,1:2),'UniformOutput',false);
    % Computing the duration in samples
    temp(:,5) = cellfun(@minus,temp(:,4),temp(:,3),'UniformOutput',false);
    % Saving blink events
    startIdx = nMsg+nSacc+1;
    endIdx = nMsg+nSacc+nBlink;
    [event(startIdx:endIdx).type] = deal('Blink');
    [event(startIdx:endIdx).sample] = temp{:,3};
    [event(startIdx:endIdx).value] = deal('');
    [event(startIdx:endIdx).duration] = temp{:,5};
    [event(startIdx:endIdx).offset] = deal(0);
    [event(startIdx:endIdx).timestamp] = temp{:,1};
end

% Extracting EFIX details
% EFIX <eye> <stime> <etime> <dur> <axp> <ayp> <aps>
if nFix > 0
    temp = regexp(fixx,'EFIX\s+[LR]\s+(\d+)\s+(\d+)\s+[0-9.]+\s+([0-9.-]+\s+[0-9.-]+)\s+[0-9.]+','tokens');
    temp = [temp{:}]';
    temp = vertcat(temp{:});
    temp(:,1:2) = cellfun(@str2num,temp(:,1:2),'UniformOutput',false);
    % Converting timestamps to samples
    temp(:,4:5) = cellfun(tstamp2sample,temp(:,1:2),'UniformOutput',false);
    % Computing the duration in samples
    temp(:,6) = cellfun(@minus,temp(:,5),temp(:,4),'UniformOutput',false);
    % Saving fixation events
    % The average x and y positions are saved as the vent value
    startIdx = nMsg+nSacc+nBlink+1;
    endIdx = nMsg+nSacc+nBlink+nFix;
    [event(startIdx:endIdx).type] = deal('Fixation');
    [event(startIdx:endIdx).sample] = temp{:,4};
    [event(startIdx:endIdx).value] = temp{:,3};
    [event(startIdx:endIdx).duration] = temp{:,6};
    [event(startIdx:endIdx).offset] = deal(0);
    [event(startIdx:endIdx).timestamp] = temp{:,1};
end

% INPUT events
if nInp > 0
    temp = {};
    temp(:,1) = cellfun(@num2str,num2cell([inp.value]'),'UniformOutput',false);
    temp(:,2) = num2cell([inp.timestamp]');
    temp(:,3) = num2cell(tstamp2sample([inp.timestamp]'));
    % Saving input events
    startIdx = nMsg+nSacc+nBlink+nFix+1;
    endIdx = nMsg+nSacc+nBlink+nFix+nInp;
    [event(startIdx:endIdx).type] = deal('Input');
    [event(startIdx:endIdx).sample] = temp{:,3};
    [event(startIdx:endIdx).value] = temp{:,1};
    [event(startIdx:endIdx).duration] = deal(1);
    [event(startIdx:endIdx).offset] = deal(0);
    [event(startIdx:endIdx).timestamp] = temp{:,2};
end

% Sorting by samples
[~,i] = sort([event.sample]);
event = event(i);

end

function asc = read_eyelink_asc_ext(filename)

% READ_EYELINK_ASC_EXT reads the header information, input triggers, messages
% and all data points from an Eyelink *.asc file
% 
% Extended by Mate Aller 21.01.2016: now eblink and sblink events are also
%   parsed. 
% 
% Use as
%   asc = read_eyelink_asc_ext(filename)

% Copyright (C) 2010-2015, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

asc.header  = {};
asc.msg     = {};
asc.input   = [];
asc.sfix    = {};
asc.efix    = {};
asc.ssacc   = {};
asc.esacc   = {};
asc.sblink  = {};
asc.eblink  = {};
asc.dat     = [];
current   = 0;

% read the whole file at once
fid = fopen(filename, 'rt');
aline = fread(fid, inf, 'char=>char');          % returns a single long string
fclose(fid);

aline(aline==uint8(sprintf('\r'))) = [];        % remove cariage return
aline = tokenize(aline, uint8(sprintf('\n')));  % split on newline

for i=1:numel(aline)
  tline = aline{i};

  if numel(tline) && any(tline(1)=='0':'9')
  % if regexp(tline, '^[0-9]')
    tmp     = sscanf(tline, '%f');
    nchan   = numel(tmp);
    current = current + 1;

    if size(asc.dat, 1)<nchan
      % increase the allocated number of channels
      asc.dat(nchan,:) = 0;
    end

    if size(asc.dat, 2)<current
      % increase the allocated number of samples
      asc.dat(:,end+10000) = 0;
    end

    % add the current sample to the data matrix
    asc.dat(1:nchan, current) = tmp;


  elseif regexp(tline, '^INPUT')
    [val, num] = sscanf(tline, 'INPUT %d %d');
    this.timestamp = val(1);
    this.value     = val(2);
    if isempty(asc.input)
      asc.input = this;
    else
      asc.input = cat(1, asc.input, this);
    end


  elseif regexp(tline, '\*\*.*')
    asc.header = cat(1, asc.header, {tline});


  elseif regexp(tline, '^MSG')
    asc.msg = cat(1, asc.msg, {tline});


  elseif regexp(tline, '^SFIX')
    asc.sfix = cat(1, asc.sfix, {tline});


  elseif regexp(tline, '^EFIX')
    asc.efix = cat(1, asc.efix, {tline});


  elseif regexp(tline, '^SSACC')
    asc.ssacc = cat(1, asc.ssacc, {tline});


  elseif regexp(tline, '^ESACC')
    asc.esacc = cat(1, asc.esacc, {tline});
    
  elseif regexp(tline, '^SBLINK')
    asc.sblink = cat(1, asc.sblink, {tline});


  elseif regexp(tline, '^EBLINK')
    asc.eblink = cat(1, asc.eblink, {tline});

  else
    % all other lines are not parsed
  end

end

% remove the samples that were not filled with real data
asc.dat = asc.dat(:,1:current);

end