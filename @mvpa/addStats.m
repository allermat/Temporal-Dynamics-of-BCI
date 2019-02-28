function mvparesObj = addStats(mvparesObj,var,varargin)
% Method for computing various statistics
%
% USAGE:
%   mvparesObj = addStats(mvparesObj,var)
%   mvparesObj = addStats(mvparesObj,var,'Name',Value)
% INPUT:
%   Required:
%       mvparesObj (object): mvpares object
%       var (string): the variable for which the statistics are computed.
%           Possible values: 'avWeights'
%   'Name'-Value arguments:
%       pathIndivFiles (cell array): full path of individual files for
%           statistics on group level data (must be specified for group 
%           level data)
%       timeWin (numeric vector): time window within the statistical 
%           analysis is performed: [starttime, endtime]. 
% OUTPUT:
%   mvparesObj (object): mvpares object with the required statistics added. 
 
% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
validVars = {'genPerf','avModelCorr','avWeights'};
addRequired(p,'mvparesObj',@(x) isa(x,'mvpares') && x.isvalid);
addRequired(p,'var',@(x) any(validatestring(x,validVars)));
addParameter(p,'pathIndivFiles',{},@(x) iscellstr(x));
addParameter(p,'timeWin',[],@(x) validateattributes(x,{'numeric'},{'vector',...
    'increasing','numel',2}));
parse(p,mvparesObj,var,varargin{:});
mvparesObj = p.Results.mvparesObj;
var = p.Results.var;
pathIndivFiles = p.Results.pathIndivFiles;
timeWin = p.Results.timeWin;

% Assign empty array and return if the result dataset is not
% generalized.
if strcmp(mvparesObj.level,'subj')
    warning('mvpa:addStats:datasetLevelMismatch',...
        ['The dataset''s level is ''subj''. Subject level statistics are ',...
        'not yet implemented. Returning.']);
    return;
else
    if isempty(pathIndivFiles)
        warning('mvpa:addStats:missingInput',...
            'pathIndivFiles must be specified for group level statistics. Returning.');
        return;
    else
        if any(~cellfun(@exist,pathIndivFiles))
            warning('mvpa:addStats:invalidInput',...
                'Each path in pathIndivFiles must point to an existing file. Returning.');
            return;
        end
    end
end

info = mvparesObj.getInfo;
% This makes sure that both Windows and Unix paths are recognized,
% regardless of the running operation system. 
temp = regexp(info.sourceFiles,'.*[\\|/](.*)\.mat','tokens','once');
sourceFileNames = cat(1,temp{:});
temp = regexp(pathIndivFiles,'.*[\\|/](.*)\.mat','tokens','once');
indivFileNames = cat(1,temp{:});

% Checking if all required datasets are specified
if any(~ismember(sourceFileNames,indivFileNames))
    warning('mvpa:addStats:missingDataset',...
        'At least one of the required individual datasets is not found. Returning. ');
    return;
end

% Loading individual objects
objList = cellfun(@mvpares,pathIndivFiles,'UniformOutput',false);

switch var
    case 'genPerf'
        fieldName = 'stat_perfEstimates';
        
        if isvector(mvparesObj.getGenTimePoints)
            genTime = 'tr';
        else
            genTime = 'tr_x_tr';
        end
        indivData = cellfun(@getGenPerfEstimates,objList,...
                            repmat({'genTime'},size(objList)),...
                            repmat({genTime},size(objList)));
        stats = mvpa.compGenPerfStats(indivData,mvparesObj.getTrTimePoints,...
                                      'genTime',genTime);

        % Setting the dataset object to writable if it is not
        if ~mvparesObj.writable, mvparesObj.setWritable(true); end
        if ismember(fieldName,fieldnames(mvparesObj.data))
            warning('mvpa:addStats:overwriteField',...
                    ['The field %s already exists in the mvpa result dataset, '...
                     'it will be overwritten.'],fieldName);
        end
        mvparesObj.data.(fieldName) = orderfields(stats);
        mvparesObj.setWritable(false);
        
    case 'avWeights'
        % Computing stats for smoothed and non smoothed weights
        for i = 1:2
            if i == 1
                smooth = false;
                fieldName = 'stat_AVmodelWeights';
            else
                smooth = true;
                fieldName = 'stat_AVmodelWeights_s';
            end
            
            indivData = cellfun(@getAVmodelWeights,objList,...
                repmat({'genTime'},size(objList)),repmat({'tr'},size(objList)),...
                repmat({'smooth'},size(objList)),repmat({smooth},size(objList)));
            if isempty(timeWin)
                stats = mvpa.compAVmodelWeightStats(indivData,...
                    mvparesObj.getTrTimePoints);
            else
                stats = mvpa.compAVmodelWeightStats(indivData,...
                    mvparesObj.getTrTimePoints,'timeWin',timeWin);
            end
            
            % Setting the dataset object to writable if it is not
            if ~mvparesObj.writable, mvparesObj.setWritable(true); end
            if ismember(fieldName,fieldnames(mvparesObj.data))
                warning('mvpa:addStats:overwriteField',...
                        ['The field %s already exists in the mvpa result dataset, '...
                         'it will be overwritten.'],fieldName);
            end
            mvparesObj.data.(fieldName) = orderfields(stats);
            mvparesObj.setWritable(false);
            
        end
    case 'avModelCorr'
        % Computing stats for smoothed and non smoothed weights
        for i = 1:2
            if i == 1
                smooth = false;
                fieldName = 'stat_AVmodelCorrelations';
            else
                smooth = true;
                fieldName = 'stat_AVmodelCorrelations_s';
            end
            genTime = 'tr';
            indivData = cellfun(@getAVmodelCorrelations,objList,...
                repmat({'all'},size(objList)),...
                repmat({'genTime'},size(objList)),repmat({genTime},size(objList)),...
                repmat({'fisherTransform'},size(objList)),repmat({true},size(objList)),...
                repmat({'smooth'},size(objList)),repmat({smooth},size(objList)));
            if isempty(timeWin)
                stats = mvpa.compAVmodelCorrStats(indivData,...
                    mvparesObj.getTrTimePoints,genTime);
            else
                stats = mvpa.compAVmodelCorrStats(indivData,...
                    mvparesObj.getTrTimePoints,genTime,'timeWin',timeWin);
            end
            % Setting the dataset object to writable if it is not
            if ~mvparesObj.writable, mvparesObj.setWritable(true); end
            if ismember(fieldName,fieldnames(mvparesObj.data))
                warning('mvpa:addStats:overwriteField',...
                        ['The field %s already exists in the mvpa result dataset, '...
                         'it will be overwritten.'],fieldName);
            end
            mvparesObj.data.(fieldName) = orderfields(stats);
            mvparesObj.setWritable(false);
            
        end
end

end