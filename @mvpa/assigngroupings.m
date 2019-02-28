function groupings = assigngroupings(cvScheme,examples_info,misc)
% Groups the examples for cross-validation

nExamples = size(examples_info,1);

if strcmp(cvScheme,'kf')
    
    % Checking input
    fieldsReqMisc = {'k','nKFrep','nTimePoints'};
    if any(~isfield(misc,fieldsReqMisc))
        error('Missing field in misc in assigngroupings! ')
    end
        
    k = misc.k;
    nKFrep = misc.nKFrep;
    nTimePoints = misc.nTimePoints;
    
    % Checking if k is smaller than the number of examples
    if k > nExamples
        error('k is greater than the number of training examples!');
    end
    
    % Generating a different K-fold grouping for each repetition of the
    % K-fold cross-validation. The groupings remain the same for across 
    % timepoints. The resulting grouping matrix has size of M x N x O 
    % with M = number of examples,  N = number of K-fold CV repetitions, 
    % O = number of time samples. 
    % One column contains numbers 1:k corresponding to the k folds. A set 
    % of examples with the same numbers corresponds to the test set for the
    % given fold. The rest of the examples correspond to the training set 
    % of the given fold. I'm also making sure, that each conditon is
    % equally well represented in each fold. 
    % Pre-allocating array for groupings
    groupings = zeros(nExamples,nKFrep);
    conditions = unique(examples_info.condition);
    if any(isnan(conditions))
        warning('Condition contains NaNs, (roughly) equal number of examples per condition per fold is not guaranteed! ');
        for i = 1:nKFrep
            groupings(:,i) = mod(randperm(nExamples)',k);
        end
    else
        % Saving the original order of the examples
        examples_info.origOrder = (1:size(examples_info,1))';
        % Function for assigning labels 0:k-1 to the input (since mod assigns
        % 0 to the k-th group).
        funRandLabel = @(x) {[x,mod(randperm(size(x,1))',k)]};
        for i = 1:nKFrep
            % Applying funRandLabel separately to the examples grouped
            % by condition and hand ensures, that each condition will
            % be roughly equally represented in each fold.
            temp = varfun(funRandLabel,examples_info,'InputVariables',...
                'origOrder','GroupingVariables',{'condition'});
            temp = cell2mat(temp.Fun_origOrder);
            temp = sortrows(temp,1);
            groupings(:,i) = temp(:,end);
        end
    end
    groupings = repmat(groupings,1,1,nTimePoints);
    % For clairaty, I change the 0 values in the grouping matrix to k. 
    groupings(groupings == 0) = k;
    
elseif strcmp(cvScheme,'loso')
    
    % Checking input
    fieldsReqMisc = {'nTimePoints'};
    if any(~isfield(misc,fieldsReqMisc))
        error('Missing field in misc in assigngroupings! ')
    end
    
    nTimePoints = misc.nTimePoints;
    % Number of predicted models per sample is the number of sessions
    % present in the dataset
    sessions = unique(examples_info.session);
    % Generating the grouping for the cross-validation resulting in a 
    % grouping matrix of M x N with M = number of examples, N = number of 
    % time points. All timepoints have the same grouping in this case. 
    groupings = examples_info.session;
    for i = 1:numel(sessions)
        groupings(examples_info.session == sessions(i)) = i;
    end
    groupings = repmat(groupings,1,nTimePoints);
     
else
    error('This cross-validation scheme is not yet implemented! ');
end

end