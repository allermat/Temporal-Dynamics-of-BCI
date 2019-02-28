function mvparesObj = mergeMvpaRes(I)
% Method for merging mvpares objects
% 
% USAGE:
%   mvparesObj = mergeMvpaRes(I)
% INPUT:
%   I (struct): structure of inputs. Required fields:
%       pathMergedFile: full path for saving the averaged file
%       pathFilesToMerge: cell array of full paths to the to be averaged
%           files
% OUTPUT:
%   mvparesObj (object): merged mvpares object

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parsing input
p = inputParser;
wrapexist = @(x) exist(x,'file');
requiredVars = {'pathMergedFile','pathFilesToMerge'};
addParameter(p,'pathMergedFile','');
addParameter(p,'pathFilesToMerge','',@(x) all(cellfun(wrapexist,x)));

parse(p,I);

pathMergedFile = p.Results.pathMergedFile;
pathFilesToMerge = p.Results.pathFilesToMerge;

if any(ismember(requiredVars,p.UsingDefaults))
    error('All required parameters must be specified!');
end

% Loading objects from disk
objList = cell(size(pathFilesToMerge));
for i = 1:numel(pathFilesToMerge)
    objList{i} = mvpares(pathFilesToMerge{i});
end
% Checking objects
checkObj = @(x) isa(x,'mvpares') && x.isvalid;
if any(~cellfun(checkObj,objList))
    error('mvpa:averageMvpaRes:invalidInput',...
          'All inputs should be valid instances of mvpares class.');
end
if numel(objList) < 2
    error('mvpa:averageMvpaRes:singleInput',...
          'There is just one input object, can''t merge, returning.');
end
if any(~ismember(cellfun(@(x) x.level,objList,'UniformOutput', ...
                         false),'subj'))
    error('mvpa:averageMvpaRes:invalidInput',...
          'All inputs should be subject level.');
end
if any(~ismember(cellfun(@(x) x.state,objList,'UniformOutput', ...
                         false),'trained_and_generalized'))
    error('mvpa:averageMvpaRes:invalidInput',...
          'All inputs should be in trained_and_generalized state.');
end
% Checking infos and detirmining merging mode (merging across
% examples is the only supported mode so far
mergeMode = checkInfos(objList);

% Merging the datasets
mergedStruct = mergeIndividualData(objList,mergeMode);

% Adding info field to the averaged dataset
indivInfos = cellfun(@getInfo,objList);
mergedInfo = objList{1}.getInfo;
infoFieldsToKeep = {'cv_scheme','gen_cond','gen_data_file','gen_data_fs',...
                    'gen_label','gen_time','sc_method','subID','svm_type','tr_cond',...
                    'tr_data_file','tr_data_fs','tr_label','tr_method','tr_timePoints'};
infoFields = fieldnames(mergedInfo);
infoFieldsToRemove = infoFields(~ismember(infoFields,infoFieldsToKeep));
mergedInfo = rmfield(mergedInfo,infoFieldsToRemove);
% The training and generalization conditions and labels are
% concatenated into one single string if necessary
fieldsToConCat = {'tr_cond','gen_cond','tr_label','gen_label'};
for i = 1:numel(fieldsToConCat)
    temp = {indivInfos.(fieldsToConCat{i})};
    if numel(unique(temp)) > 1
        mergedInfo.(fieldsToConCat{i}) = strjoin(temp,'_&_');
    else
        mergedInfo.(fieldsToConCat{i}) = temp{1};
    end
end
% The training and generalization examples are the union of the
% corresponding individual examples
temp = {indivInfos.tr_isExample};
temp = cat(2,temp{:});
mergedInfo.tr_isExample = any(temp,2);
temp = {indivInfos.gen_isExample};
temp = cat(2,temp{:});
mergedInfo.gen_isExample = any(temp,2);
% cv_scheme specific fields
if strcmp(mergedInfo.cv_scheme,'kf')
    mergedInfo.k = sum([indivInfos.k]);
    % nKFrep should be uniform across datasets
    mergedInfo.nKFrep = indivInfos(1).nKFrep;
end
mergedInfo.sourceFiles = pathFilesToMerge;
mergedInfo = orderfields(mergedInfo);
mergedStruct.info = mergedInfo;
mergedStruct = orderfields(mergedStruct); %#ok<NASGU>
% Saving merged dataset
if exist(pathMergedFile,'file')
    warning('The file to be written already exists, overwriting. ');
end

save(pathMergedFile,'-struct','mergedStruct','-v7.3');
% Loading mvpares object from the saved dataset
mvparesObj = mvpares(pathMergedFile);

end

function mergeMode = checkInfos(objList)

infos = cellfun(@getInfo,objList);

if numel(unique({infos.cv_scheme})) > 1
    error('mvpa:mergeMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The cv_scheme is not consistent across mvpa results for merging.');
elseif numel(unique({infos.subID})) > 1
    error('mvpa:mergeMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The subID is not consistent across mvpa results for merging.');
elseif numel(unique({infos.gen_time})) > 1
    error('mvpa:mergeMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The gen_time is not consistent across mvpa results for merging.');
elseif numel(unique({infos.sc_method})) > 1
    error('mvpa:mergeMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The sc_method is not consistent across mvpa results for merging.');
elseif numel(unique({infos.svm_type})) > 1
    error('mvpa:mergeMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The svm_type is not consistent across mvpa results for merging.');
elseif numel(unique({infos.tr_method})) > 1
    error('mvpa:mergeMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_method is not consistent across mvpa results for merging.');
elseif size(unique(cell2mat(cat(1,infos.tr_timePoints)),'rows'),1) > 1
    error('mvpa:mergeMvpaRes:checkInfos:inconsistentMvpaResults',...
        'The tr_timePoints are not consistent across mvpa results for merging.');
end

if isfield(infos,'gen_data_file')
    if numel(unique({infos.gen_data_file})) > 1
        error('mvpa:mergeMvpaRes:checkInfos:inconsistentMvpaResults',...
              'The gen_data_file is not consistent across mvpa results for merging.');
    end
end

if all(ismember({infos.cv_scheme},'kf'))
    if numel(unique([infos.nKFrep])) > 1
        error('mvpa:mergeMvpaRes:checkInfos:inconsistentMvpaResults',...
              'The nKFrep is not consistent across mvpa results for merging.');
    end
end

% Checking if there are different generalization conditions
if numel(unique({infos.gen_cond})) ~= numel({infos.gen_cond})
    warning('mvpa:mergeMvpaRes:checkInfos:dependentGenExamples',...
        ['Not all of the mvpa results'' gen_cond is distinct, some ' ...
         'of the resulting generalized examples will be dependent.']);
end
mergeMode = 'examples';

end

function mergedStruct = mergeIndividualData(objList,mergeMode)
% Function for merging individual level data

if strcmp(mergeMode,'examples')
    
    mergedFieldNames = {'gen_accuracies','gen_examples','gen_groupings', ...
                        'gen_predlabels','tr_examples','tr_accuracies'};
    dimToMerge = [2,1,1,2,1,2];
    
    mergedStruct = struct();
    
    for i = 1:numel(mergedFieldNames)
        temp = cellfun(@(x) x.data.(mergedFieldNames{i}),objList, ...
                'UniformOutput',false);
        
        if strcmp(mergedFieldNames{i},'gen_groupings')
            % The total number of folds in the merged dataset will be
            % the sum of the folds over the individual datasets
            for j = 1:numel(temp)
                if j == 1, continue; end
               temp{j} = temp{j} + max(temp{j-1}(:));
            end
        elseif ismember(mergedFieldNames{i},{'gen_examples','tr_examples'})
            % We have to take care of the possibility, that the
            % variables in the example info tables are not
            % consistent across datasets
            allFieldNames = cellfun(@(x) ...
                                    x.Properties.VariableNames,temp,'UniformOutput',false);
            allFieldNames = unique([allFieldNames{:}]);
            for j = 1:numel(temp)
                % Adding vectors of NaNs to the example info
                % tables in place of the missing variables
                id = ~ismember(allFieldNames, ...
                               temp{j}.Properties.VariableNames);
                missingVarNames = allFieldNames(id);
                t = NaN(size(temp{j},1),numel(missingVarNames));
                t = cell2table(num2cell(t),'VariableNames',missingVarNames);
                temp{j} = cat(2,temp{j},t);
                % Adding an extra variable to the example info tables
                % to keep track which examples came from which dataset
                temp{j}.mergedDataSetId = ones(size(temp{j},1),1)*j;
            end
        end
        % Concatenating the required fields along the specified dimensions
        if isa(temp{1},'double')
            mergedStruct.(mergedFieldNames{i}) = ...
                catpad(dimToMerge(i),temp{:});
        else
            mergedStruct.(mergedFieldNames{i}) = ...
                cat(dimToMerge(i),temp{:});
        end
    end
    
end

end