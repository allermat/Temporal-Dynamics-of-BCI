function summary = permttest(data,H0,varargin)
% One sample bootstrap test based on t statistic 
% 
% SYNTAX: 
%         summary = permttest(data,H0)
%         summary = permttest(data,H0,'Name','Value')
% DETAILS: 
%          By default it does a two-sided hypothesis test. 
%          NaN safe. 
%          cf. Efron, Tibshirani introduction to bootstrap. 
% INPUT: 
%       data: vector of sample of data points
%       H0: the value against witch the one sample T-test is to be
%           performed (eg. the null hypothesis);
%       Name-value pair arguments:
%           nPerm: number of permutations to be drawn. 
%                  Default = 5000. 
%           tail: type of alternative hypothesis to evaluate. 
%                 'both': Test the alternative hypothesis that the 
%                         population mean is not m.
%                 'left': Test the alternative hypothesis that the 
%                         population mean is less than m.
%                 'right': Test the alternative hypothesis that the 
%                          population mean is greater than m. 
%                 Default = 'both'. 
%           'useParallel': If true and if a parpool of the 
%               Parallel Computing Toolbox™ is open, compute bootstrap 
%               iterations in parallel. If the Parallel Computing Toolbox 
%               is not installed, or a parpool is not open, computation 
%               occurs in serial mode. Default is false, meaning serial 
%               computation.
% OUTPUT: 
%       Dataset array with following fields: 
%       p: the probability of observing a test statistic as extreme as, 
%          or more extreme than, the observed value under the null 
%          hypothesis.
%       t: the t obtained from the test. 
%       df: degrees of freedom. 
%

%% Parsing input parameters. 
parser = inputParser;

checknBootstrap = @(x) isscalar(x) && x > 0;
validTails = {'both','left','right'};
checkTail = @(x) any(validatestring(x,validTails));

addRequired(parser,'data',@iscolumn);
addRequired(parser,'H0',@isscalar);
addParameter(parser,'nPerm',5000,checknBootstrap);
addParameter(parser,'tail','both',checkTail);
addParameter(parser,'useParallel',false,@islogical);

parse(parser,data,H0,varargin{:});

nPerm = parser.Results.nPerm;
tail = parser.Results.tail;
useParallel = parser.Results.useParallel;

%% 
% Removing NaNs and Infs. 
data = data(isfinite(data));

if length(data) > 1  % We need at least 2 data points to do bootstrapping. 
    % One sample ttest against mu.
    funTtest = @(x,mu) sqrt(length(x))*(mean(x)-mu)/std(x);
    
    tOrig = funTtest(data,H0);
    df = length(data)-1;
    
    indivPermT = NaN(nPerm,1);
    
    if useParallel
        parfor i = 1:nPerm
            % Generating random numbers
            perm = randn(size(data));
            % Dividing them by their absolute values to get +/- 1s.
            perm = perm./abs(perm);
            % Multiplying the original biases with the permutation vector, element
            % by element.
            indivPermT(i) = funTtest(data.*perm,mean(data));
        end
    else
        for i = 1:nPerm
            % Generating random numbers
            perm = randn(size(data));
            % Dividing them by their absolute values to get +/- 1s.
            perm = perm./abs(perm);
            % Multiplying the original biases with the permutation vector, element
            % by element.
            indivPermT(i) = funTtest(data.*perm,mean(data));
        end
    end
    
    switch tail
        case 'both'
            p = 2*min(sum(indivPermT <= tOrig),...
                 sum(indivPermT >= tOrig))/nPerm;
        case 'left'
            p = sum(indivPermT <= tOrig)/nPerm;
        case 'right'
            p = sum(indivPermT >= tOrig)/nPerm;
    end
    
else
    tOrig = NaN;
    df = NaN;
    p = NaN;
end

summary = table(p,tOrig,df,'VariableNames',{'p','t','df'});

end