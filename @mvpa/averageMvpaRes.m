function mvparesObj = averageMvpaRes(I)
% Method for averaging mvpares objects
% 
% USAGE:
%   mvparesObj = averageMvpaRes(I)
% INPUT:
%   I (struct): structure of inputs. Required fields:
%       pathAveragedFile: full path for saving the averaged file
%       pathFilesToAverage: cell array of full paths to the to be averaged
%           files
% OUTPUT:
%   mvparesObj (object): averaged mvpares object

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
wrapexist = @(x) exist(x,'file');
requiredVars = {'pathAveragedFile','pathFilesToAverage'};
addParameter(p,'pathAveragedFile','');
addParameter(p,'pathFilesToAverage','',@(x) all(cellfun(wrapexist,x)));

parse(p,I);

pathAveragedFile = p.Results.pathAveragedFile;
pathFilesToAverage = p.Results.pathFilesToAverage;

if any(ismember(requiredVars,p.UsingDefaults))
    error('All required parameters must be specified!');
end

% Loading objects from disk
objList = cell(size(pathFilesToAverage));
for i = 1:numel(pathFilesToAverage)
    objList{i} = mvpares(pathFilesToAverage{i});
end
% Checking objects
checkObj = @(x) isa(x,'mvpares') && x.isvalid;
if any(~cellfun(checkObj,objList))
    error('mvpa:averageMvpaRes:invalidInput',...
        'All inputs should be valid instances of mvpares class.');
end
if numel(objList) < 2
    error('mvpa:averageMvpaRes:singleInput',...
        'There is just one input object, can''t average, returning.');
end
checkInfos(objList);

% Looking for data available for averaging in the objects
dataFields = cellfun(@who,objList,'UniformOutput',false);
averagedDataNames = {'gen_AVmodelCorrelations','gen_AVmodelCorrelations_s',...
    'gen_AVmodelEstimates','gen_AVmodelEstimates_s','gen_AVmodelWeights',...
    'gen_AVmodelWeights_s','gen_perfEstimates'};
averages = cell(size(averagedDataNames));

for i = 1:numel(averagedDataNames)
    
    switch averagedDataNames{i}
        case 'gen_AVmodelCorrelations'
            indivDataName = 'gen_AVmodelCorrelations';
        case 'gen_AVmodelCorrelations_s'
            indivDataName = 'gen_AVmodelCorrelations';
        case 'gen_AVmodelEstimates'
            indivDataName = 'gen_AVmodelEstimates';
        case 'gen_AVmodelEstimates_s'
            indivDataName = 'gen_AVmodelEstimates';
        case 'gen_AVmodelWeights'
            indivDataName = 'gen_AVmodelEstimates';
        case 'gen_AVmodelWeights_s'
            indivDataName = 'gen_AVmodelEstimates';
        case 'gen_perfEstimates'
            indivDataName = 'gen_perfEstimates';
    end
    
    if all(cellfun(@ismember,repmat({indivDataName},size(dataFields)),dataFields))
        if any(~cellfun(@ismember,repmat({indivDataName},size(dataFields)),dataFields))
            warning('mvpa:averageMvpaRes:missingData',...
                ['Not all of the input objects have the ''%s'' ',...
                'field. This field will not be averaged.'],indivDataName);
        else
            averages{i} = averageIndividualData(objList,averagedDataNames{i});
        end
    end
end

if any(~cellfun(@isempty,averages))
    keepIdx = ~cellfun(@isempty,averages);
    % Removing fields which are not averaged
    averages = averages(keepIdx);
    fieldsAveraged = averagedDataNames(keepIdx);
    avgStruct = cell2struct(averages,fieldsAveraged,2);
    % Adding info field to the averaged dataset
    avgInfo = objList{1}.getInfo;
    infoFieldsToKeep = {'cv_scheme','gen_cond','gen_data_file','gen_data_fs',...
        'gen_label','gen_time','sc_method','subID','svm_type','tr_cond',...
        'tr_data_file','tr_data_fs','tr_label','tr_method','tr_timePoints'};
    infoFields = fieldnames(avgInfo);
    infoFieldsToRemove = infoFields(~ismember(infoFields,infoFieldsToKeep));
    avgInfo = rmfield(avgInfo,infoFieldsToRemove);
    avgInfo.subID = 'group';
    avgInfo.sourceFiles = pathFilesToAverage;
    avgInfo = orderfields(avgInfo);
    avgStruct.info = avgInfo;
    avgStruct = orderfields(avgStruct); %#ok<NASGU>
    % Saving averaged dataset
    save(pathAveragedFile,'-struct','avgStruct','-v7.3');
    % Loading mvpares object from the saved dataset
    mvparesObj = mvpares(pathAveragedFile);
else
    warning('mvpa:averageMvpaRes:requestedDataNotAvailable',...
        'Could not find any data to average. Retruning ');
    mvparesObj = [];
end

end

function checkInfos(objList)

infos = cellfun(@getInfo,objList);

if numel(unique({infos.cv_scheme})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The cv_scheme is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.gen_cond})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The gen_cond is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.gen_label})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The gen_label is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.gen_time})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The gen_time is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.sc_method})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The sc_method is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.svm_type})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The svm_type is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.tr_cond})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_cond is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.tr_label})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_label is not consistent across mvpa results for averaging.');
elseif numel(unique({infos.tr_method})) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_method is not consistent across mvpa results for averaging.');
elseif size(unique(cell2mat(cat(1,infos.tr_timePoints)),'rows'),1) > 1
    error('mvpa:averageMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_timePoints are not consistent across mvpa results for averaging.');
end

end

function avgData = averageIndividualData(objList,dataType)
% Function for averaging individual level data

temp = objList{1}.getInfo;
genTime = temp.gen_time;

switch dataType
    case 'gen_AVmodelCorrelations'
        indivData = cellfun(@getAVmodelCorrelations,objList,...
            repmat({'all'},size(objList)),...
            repmat({'genTime'},size(objList)),repmat({genTime},size(objList)));
        fieldExclusionStr = '';
    case 'gen_AVmodelCorrelations_s'
        indivData = cellfun(@getAVmodelCorrelations,objList,...
            repmat({'all'},size(objList)),...
            repmat({'genTime'},size(objList)),repmat({'tr'},size(objList)),...
            repmat({'smooth'},size(objList)),repmat({true},size(objList)));
        fieldExclusionStr = '';
    case 'gen_AVmodelEstimates'
        indivData = cellfun(@getAVmodelEstimates,objList,...
            repmat({'genTime'},size(objList)),repmat({genTime},size(objList)));
        fieldExclusionStr = '_ci$';
    case 'gen_AVmodelEstimates_s'
        indivData = cellfun(@getAVmodelEstimates,objList,...
            repmat({'genTime'},size(objList)),repmat({'tr'},size(objList)),...
            repmat({'smooth'},size(objList)),repmat({true},size(objList)));
        fieldExclusionStr = '_ci$';
    case 'gen_AVmodelWeights'
        indivData = cellfun(@getAVmodelWeights,objList,...
            repmat({'genTime'},size(objList)),repmat({genTime},size(objList)));
        fieldExclusionStr = '_ci$';
    case 'gen_AVmodelWeights_s'
        indivData = cellfun(@getAVmodelWeights,objList,...
        repmat({'genTime'},size(objList)),repmat({'tr'},size(objList)),...
        repmat({'smooth'},size(objList)),repmat({true},size(objList)));
        fieldExclusionStr = '_ci$';
    case 'gen_perfEstimates'
        indivData = cellfun(@getGenPerfEstimates,objList,...
            repmat({'genTime'},size(objList)),repmat({genTime},size(objList)));
        fieldExclusionStr = '_perFolds$';
end

avgData = indivData(1);
fieldNames = fieldnames(avgData);
fieldsToAverage = fieldNames(cellfun(@isempty,regexp(fieldNames,fieldExclusionStr)));
avgData = rmfield(avgData,fieldNames(...
    ~cellfun(@isempty,regexp(fieldNames,fieldExclusionStr))));

for i = 1:numel(fieldsToAverage)
    dimToCat = sum(size(indivData(1).(fieldsToAverage{i})) > 1)+1;
    temp = cat(dimToCat,indivData.(fieldsToAverage{i}));
    
    if ~isempty(regexp(fieldsToAverage{i},'_beta$','once'))
        avgData.(fieldsToAverage{i}) = mean(temp,dimToCat);
        % Adding group level error estimates
        fieldName = strrep(fieldsToAverage{i},'_beta','_std');
        avgData.(fieldName) = std(temp,0,dimToCat);
        fieldName = strrep(fieldsToAverage{i},'_beta','_err');
        avgData.(fieldName) = (std(temp,0,dimToCat))./sqrt(size(temp,3));
        fieldName = strrep(fieldsToAverage{i},'_beta','_ci');
        [ciLow,ciHigh] = mvpa.confmean(temp,dimToCat,0.32);
        avgData.(fieldName) = (ciHigh-ciLow)/2;
    elseif ~isempty(regexp(fieldsToAverage{i},'_wav$','once'))
        avgData.(fieldsToAverage{i}) = circ_mean(temp,[],dimToCat);
        % Adding group level ci
        fieldName = strrep(fieldsToAverage{i},'_wav','_ci');
        warning('off','all');
        avgData.(fieldName) = circ_confmean(temp,0.32,[],[],dimToCat);
        warning('on','all');
    elseif ~isempty(regexp(fieldsToAverage{i},'^r2_','once'))
        % Fisher transform first, then average when averaging correlation
        % coefficients.
        avgData.(fieldsToAverage{i}) = ...
            squeeze((tanh(mean(atanh(sqrt(temp)),dimToCat))).^2);
    elseif ~isempty(regexp(fieldsToAverage{i},'^r_','once'))
        % Fisher transform first, then average when averaging correlation
        % coefficients.
        avgData.(fieldsToAverage{i}) = ...
            squeeze(tanh(mean(atanh(temp),dimToCat)));
    elseif ~isempty(regexp(fieldsToAverage{i},'^(acc|b|int|rf)_','once'))
        avgData.(fieldsToAverage{i}) = mean(temp,dimToCat);
        % Adding group level error estimates
        fieldName = [fieldsToAverage{i},'_std'];
        avgData.(fieldName) = std(temp,0,dimToCat);
        fieldName = [fieldsToAverage{i},'_err'];
        avgData.(fieldName) = (std(temp,0,dimToCat))./ ...
            sqrt(size(temp,3));
    elseif ~isempty(regexp(fieldsToAverage{i},'^cc_','once'))
        % Fisher transform first, then average when averaging correlation
        % coefficients.
        avgData.(fieldsToAverage{i}) = ...
            squeeze(tanh(mean(atanh(temp),dimToCat)));
        % Adding group level error estimates
        fieldName = strrep(fieldsToAverage{i},'cc_','std_');
        avgData.(fieldName) = squeeze(tanh(std(atanh(temp),0,dimToCat)));
        fieldName = strrep(fieldsToAverage{i},'cc_','err_');
        avgData.(fieldName) = squeeze(tanh(std(atanh(temp),0,dimToCat)))./sqrt(size(temp,3));
        fieldName = strrep(fieldsToAverage{i},'cc_','ci_');
        [tempLow,tempHigh] = mvpa.confmean(atanh(temp),dimToCat,0.32);
        ciLow = squeeze(tanh(tempLow));
        ciHigh = squeeze(tanh(tempHigh));
        avgData.(fieldName) = (ciHigh-ciLow)/2;
    end
    
end
avgData = orderfields(avgData);

end