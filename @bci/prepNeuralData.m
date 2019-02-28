function bciData = prepNeuralData(info,actSample,varargin)
% Converts predicted labels to the format accepatble for the modelfitting

% Parsing input
p = inputParser;
validFoldActions = {'avg','none','shrink'};
addRequired(p,'info',@(x)validateattributes(x,{'table'},...
    {'nonempty'}));
addRequired(p,'actSample',@(x)validateattributes(x,{'numeric'},...
    {'size',[size(info,1),1]}));
addOptional(p,'foldAction','none',@(x) any(validatestring(x,validFoldActions)));
parse(p,info,actSample,varargin{:});
info = p.Results.info;
actSample = p.Results.actSample;
foldAction = p.Results.foldAction;

info.actSample = actSample;

switch foldAction
    case 'avg'
        % Averaging over examples of conditons within each fold
        temp = varfun(@nanmean,info,'InputVariables',{'actSample'},...
            'GroupingVariables',{'session','task','locA','locV','relV'});
        temp.Properties.VariableNames{'nanmean_actSample'} = 'actSample';
        temp.Properties.RowNames = {};
    case 'shrink'
        % Taking just those examples from each fold, which correspond to
        % the session from the cross-validation
        temp = info(info.fold == info.session,:);
    otherwise
        temp = info;
end

[temp.respV,temp.respA] = deal(NaN(size(temp,1),1));
temp.respV(temp.task == 2) = temp.actSample(temp.task == 2);
temp.respA(temp.task == 1) = temp.actSample(temp.task == 1);

% Selecting appropriate trials and variables for BCI fitting
bciData = temp(...
    ~isnan(temp.locV) & ...    % Audio-visual trials
    ~isnan(temp.locA), ...
    {'locV','locA','relV','respV','respA'});  % and these variables

end

