function labelOut = selectlabel(inputLabel,inputConds)
% Selects the label for the given input label and condition
% 
% INPUT:
%   inputLabel = input label
%   inputCond = input condition
% 
% OUTPUT:
%   labelOut = the output label (for training or generalization)
% 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

%% Parsing input
p = inputParser;

validLabels = {'disp','finger','hand','resp','stim','sAhatInd','sVhatInd','sHatComm',...
    'sAhatANDsVhat','task','relV','LvsRa','LvsRv'};
validConds = {'A','V','AV','AxorV','pseudoAV','c','ci','i','a','v','av','hrel','lrel',...
    'lh','rh','la','ha'};

addRequired(p,'inputLabel',@(x)any(validatestring(x,validLabels)));
addRequired(p,'inputConds',@(x)all(ismember(regexp(x,'-','split'),validConds)));

% Parsing inputs.
parse(p,inputLabel,inputConds);

% Assigning inputs to variables
inputLabel = p.Results.inputLabel;
inputConds = p.Results.inputConds;

%% Main
% Splitting the input conditions to parts
conds = regexp(inputConds,'-','split');

if strcmp(inputLabel,'stim')
    if ismember('a',conds)
        labelOut = 'locA';
    elseif ismember('v',conds)
        labelOut = 'locV';
    elseif ismember('av',conds)
        labelOut = 'locV';
    elseif ismember('A',conds)
        labelOut = 'locA';
    elseif ismember('V',conds)
        labelOut = 'locV';
    elseif ismember('AV',conds)
        labelOut = 'locV';
    elseif ismember('pseudoAV',conds)
        labelOut = 'locV';
    elseif ismember('AxorV',conds)
        labelOut = 'locAxorV';
    else
        error('Unidentified condition');
    end
elseif strcmp(inputLabel,'resp')
    labelOut = 'resp';
elseif strcmp(inputLabel,'relV')
    labelOut = 'relV';
elseif strcmp(inputLabel,'disp')
    labelOut = 'disp';
elseif strcmp(inputLabel,'task')
    labelOut = 'task';
elseif strcmp(inputLabel,'hand')
    labelOut = 'hand';
elseif strcmp(inputLabel,'finger')
    labelOut = 'finger';
elseif strcmp(inputLabel,'sAhatInd')
    labelOut = 'sA_hat_indep';
elseif strcmp(inputLabel,'sVhatInd')
    labelOut = 'sV_hat_indep';
elseif strcmp(inputLabel,'sHatComm')
    labelOut = 's_hat_comm';
elseif strcmp(inputLabel,'sAhatANDsVhat')
    labelOut = 'sA_hat_and_sV_hat';
elseif strcmp(inputLabel,'LvsRa')
    labelOut = 'LvsRa';
elseif strcmp(inputLabel,'LvsRv')
    labelOut = 'LvsRv';
else
    error('Unidentified label');
end

end