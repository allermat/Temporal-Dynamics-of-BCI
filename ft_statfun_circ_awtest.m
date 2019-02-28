function [s,cfg] = ft_statfun_circ_awtest(cfg,dat,design)
% Wrapper function for circ_AWTest to be used with fieldtrip statistic functions
% 
% DETAILS: 
%   Calculates the Anderson-Wu test on the biological data in dat (the 
%   dependent variable), using the information on the independent variable 
%   (ivar) in design. The function uses circ_AWTest implemented by
%   Tim Rohe based on Anderson, Wu (1995)
% USAGE:
%   Use this function by calling one of the high-level statistics functions
%       [stat] = ft_timelockstatistics(cfg, timelock1, timelock2, ...)
%       [stat] = ft_freqstatistics(cfg, freq1, freq2, ...)
%       [stat] = ft_sourcestatistics(cfg, source1, source2, ...)
%       with the following configuration option
%       cfg.statistic = 'ft_statfun_depsamplesFunivariate'
%       see FT_TIMELOCKSTATISTICS, FT_FREQSTATISTICS or FT_SOURCESTATISTICS
%       for details.
%   For low-level use, the external interface of this function has to be
%       [s,cfg] = ft_statfun_circ_awtest(cfg,dat,design);
% INPUT:
%   dat (numerical): contains the biological data, Nsamples x Nreplications
%   design (numerical): contains the independent variable (ivar) and the 
%       unit-of-observation (uvar)
%       factor, Nfac x Nreplications
%   cfg (struct): configuration structure, with the following options
%       cfg.computestat = 'yes' or 'no', calculate the statistic 
%           (default='yes')
%       cfg.computecritval = 'yes' or 'no', calculate the critical values 
%           of the test statistics (default='no')
%       cfg.computeprob  = 'yes' or 'no', calculate the p-values 
%           (default='no')
%       The following options are relevant if 
%           cfg.computecritval='yes' and/or cfg.computeprob='yes'.
%       cfg.alpha = critical alpha-level of the statistical test 
%           (default=0.05)
%       cfg.tail  = -1, 0, or 1, left, two-sided, or right (default=1)
%       cfg.tail in combination with cfg.computecritval='yes'
%           determines whether the critical value is computed at
%           quantile cfg.alpha (with cfg.tail=-1), at quantiles
%           cfg.alpha/2 and (1-cfg.alpha/2) (with cfg.tail=0), or at
%           quantile (1-cfg.alpha) (with cfg.tail=1). For the
%           Fstatistic only cfg.tail = 1 makes sense.
%       Design specification
%       cfg.ivar  = row number of the design that contains the labels of 
%           the conditions that must be  compared (default=1). The labels 
%           range from 1 to the number of conditions.
%       cfg.uvar  = row number of design that contains the labels of the 
%           units-of-observation (subjects or trials) (default=2). The 
%           labels are assumed to be integers ranging from 1 to the number 
%           of units-of-observation.
%       cfg.whicheffect = the effect of interest. Possible values: 
%           {'A', 'B', 'C', 'A X B', 'A X C', 'B X C', 'A X B X C'}
% OUTPUT:
%   s (sturcture): results of the statistical test. Contains the following
%       fields depending on input settings:
%       s.stat (vector): nSamples x 1 array of F-values
%       s.prob (vector): nSamples x 1 array of p-values
%       s.critval (scalar): Critical value of the F-statisic at the 
%           specified alpha level

% Copyright (C) 2018, Mate Aller
% allermat@gmail.com
%
% Adapted from ft_statfun_depsamplesFunivariate
% Copyright (C) 2014, Diego Lozano-Soldevilla
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
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


% set defaults
if ~isfield(cfg, 'computestat'),       cfg.computestat='yes';     end;
if ~isfield(cfg, 'computecritval'),    cfg.computecritval='no';   end;
if ~isfield(cfg, 'computeprob'),       cfg.computeprob='no';      end;
if ~isfield(cfg, 'alpha'),             cfg.alpha=0.05;            end;
if ~isfield(cfg, 'tail'),              cfg.tail=1;                end;

conds = unique(design(cfg.ivar,:)', 'rows');
nconds = length(conds);

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
    error('P-values can only be calculated if the test statistics are calculated.');
end

% perform some checks on the design
nuospercond=zeros(nconds,1);
for condindx=1:nconds
    nuospercond(condindx) = sum(ismember(design(cfg.ivar,:)',conds(condindx,:),'rows'));
end

nsmpls = size(dat,1);
nFactors = numel(cfg.ivar);

if strcmp(cfg.computestat,'yes')
    
    if ~isfield(cfg, 'whicheffect') 
        error('whicheffect must be specified!');
    end
    
    % compute the statistic
    stat = zeros(nsmpls,1);
    prob = zeros(nsmpls,1);
    whicheffect = cfg.whicheffect;
    if nFactors == 1
        grpNames = {'A'};
        grpA = design(cfg.ivar(1),:)';
        grpB = [];
        grpC = [];
    elseif nFactors == 2
        grpNames = {'A','B'};
        grpA = design(cfg.ivar(1),:)';
        grpB = design(cfg.ivar(2),:)';
        grpC = [];
    elseif nFactors == 3
        grpNames = {'A','B','C'};
        grpA = design(cfg.ivar(1),:)';
        grpB = design(cfg.ivar(2),:)';
        grpC = design(cfg.ivar(3),:)';
    else
        error('The test is implemented only for one-way(2 factors), 2x2 and 2x2x2 designs');
    end
    
    parfor smplindx = 1:nsmpls
        datonesmpl = dat(smplindx,:)';
        
        tbl = circ_AWTest(datonesmpl,grpA,grpB,grpC,grpNames);
        
        % Extracting the LRTS for the required effect from the table
        stat(smplindx) = tbl{ismember(tbl(:,1),{whicheffect}),2};
        prob(smplindx) = tbl{ismember(tbl(:,1),{whicheffect}),3};
    end
    s.stat = stat;
end

if strcmp(cfg.computecritval,'yes')
    % also compute the critical values
    % LRTS (i.e. difference between log likelihood of 2 models) is
    % approximately chi^2 distributed (Wilkin's theorem, see
    % http://en.wikipedia.org/wiki/Likelihood-ratio_test)
    % We estimate one mean for each factor for the  alternative
    % model and one grand mean for null model
    s.dfnum = nFactors - 1;
    if cfg.tail == -1
        error(['For a dependent samples LRTS-statistic, it does not ' ...
               'make sense to calculate a left tail critical value.']);
    end
    if cfg.tail == 0
        error(['For a dependent samples LRTS-statistic, it does not ' ...
               'make sense to calculate a two-sided critical value.']);
    end
    if cfg.tail == 1
        s.critval = icdf('Chisquare',cfg.alpha,s.dfnum);
    end
end

if strcmp(cfg.computeprob,'yes')
    % also return the p-values
    if cfg.tail==-1
        error(['For a dependent samples LRTS-statistic, it does not ' ...
               'make sense to calculate a left tail p-value.']);
    end
    if cfg.tail==0
        error(['For a dependent samples LRTS-statistic, it does not ' ...
               'make sense to calculate a two-sided p-value.']);
    end
    if cfg.tail==1
        s.prob = prob;
    end
end
