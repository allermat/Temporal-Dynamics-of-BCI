function [statTab] = circ_AWTest(theta,grpA,grpB,grpC,grpNames);
% computes significance of main and interaction effects (only in 2x2 and 2x2x2 designs) on circular data,
%based on a log likelihood ratio statistic. Described in Anderson, Wu
%(1995) Measuring location effects from factorial experiments with a
%directional response
% Contribution of Tim Rohe
%
% In:       theta     column vector of angles in radians
%           grpA / B / C column vector of indices indicating group
%                     membership of factors A, B, C (leave C empty for 2 way design)
%           grpNames  cell containing group names    
% out:      statTab   LRTS statistic and p values for main and interaction effects 
if isempty(grpB);
    nFact = 1;
    if length(unique(grpA))~=2 ;
        error('Only for two groups implemented!');
    end;
elseif isempty(grpC);
    nFact = 2;
    if length(unique(grpA))~=2 || length(unique(grpB))~= 2;
        error('Only for 2x2 designs');
    end;

else
    nFact = 3;
    if length(unique(grpA))~=2 || length(unique(grpB))~= 2 || length(unique(grpC))~= 2;
        error('Only for 2x2x2 designs');
    end;
end;


% compute maximum likelihood estimate of concentration parameter k as
% described by Fisher 1995 and implemented by Berens 2009
[muHat kHat] = circ_vmpar(theta);
theta_grand = computeMeanTheta_resultantLength(theta,ones(length(theta),1));
% (re-)code index variables for main and interaction effects
n = length(theta);
grpATmp = zeros(n,1);
valA = unique(grpA);
grpATmp(grpA == valA(1)) = 1;
grpATmp(grpA == valA(2)) = 2;
grpA = grpATmp;
if nFact == 2 | nFact == 3;
    grpBTmp = zeros(n,1);
    valB = unique(grpB);
    grpBTmp(grpB == valB(1)) = 1;
    grpBTmp(grpB == valB(2)) = 2;
    grpB = grpBTmp;
    % 2 way interaction effect index by surrogate variable
    grpAxB = zeros(n,1);
    ind = grpA ~= grpB;
    grpAxB(ind) = 1;
    ind = grpA == grpB;
    grpAxB(ind) = 2;
end;
if nFact == 3;
    grpCTmp = zeros(n,1);
    valC = unique(grpC);
    grpCTmp(grpC == valC(1)) = 1;
    grpCTmp(grpC == valC(2)) = 2;
    grpC = grpCTmp;
    %-- left 2 way interaction effect indices
    grpAxC = zeros(n,1);
    ind = grpA ~= grpC;
    grpAxC(ind) = 1;
    ind = grpA == grpC;
    grpAxC(ind) = 2;
    grpBxC = zeros(n,1);
    ind = grpB ~= grpC;
    grpBxC(ind) = 1;
    ind = grpB == grpC;
    grpBxC(ind) = 2;
    
    %-- 3 way interaction
    grpAxBxC = grpAxB;
    % reverse interaction for C == 1 vs. C == 2, so that one group has odd
    % and other has even number of high levels (cf. p. 14)
    grpAxBxC((grpC == 2) & (grpAxB == 1)) = 2;
    grpAxBxC((grpC == 2) & (grpAxB == 2)) = 1;    
end;

%-- compute log likelihood ratio test statistic-----------
%-- Main effect A
[theta_A,r_A] = computeMeanTheta_resultantLength(theta,grpA);
LRTS_A = 2*kHat * sum(r_A .* (1 - cos(theta_A - theta_grand)));

%--- compute effects size (Nagelkerkes coefficient of determination) from
%LRTS ( compare Nagelkerkes equation 1b and Anderson&Wu s appendix 1)
R2_A = 1-exp(-LRTS_A / n);

if nFact == 2 | nFact == 3;
    %-- Main effect B
    [theta_B,r_B] = computeMeanTheta_resultantLength(theta,grpB);
    LRTS_B = 2*kHat * sum(r_B .* (1 - cos(theta_B - theta_grand)));
    R2_B = 1-exp(-LRTS_B / n);
    
    %-- Interaction AxB
    [theta_AxB,r_AxB] = computeMeanTheta_resultantLength(theta,grpAxB);
    LRTS_AxB = kHat * sum(r_AxB .* (1 - abs(cos(theta_AxB - theta_grand))));
    R2_AxB = 1-exp(-LRTS_AxB / n);
end;
if nFact == 3;
    %-- Main effect C
    [theta_C,r_C] = computeMeanTheta_resultantLength(theta,grpC);
    LRTS_C = 2*kHat * sum(r_C .* (1 - cos(theta_C - theta_grand)));
    R2_C = 1-exp(-LRTS_C / n);
    
    %-- Interaction AxC
    [theta_AxC,r_AxC] = computeMeanTheta_resultantLength(theta,grpAxC);
    LRTS_AxC = kHat * sum(r_AxC .* (1 - abs(cos(theta_AxC - theta_grand))));
    R2_AxC = 1-exp(-LRTS_AxC / n);
    
    %-- Interaction BxC
    [theta_BxC,r_BxC] = computeMeanTheta_resultantLength(theta,grpBxC);
    LRTS_BxC = kHat * sum(r_BxC .* (1 - abs(cos(theta_BxC - theta_grand))));
    R2_BxC = 1-exp(-LRTS_BxC / n);
    
    %-- Interaction AxBxC
    [theta_AxBxC,r_AxBxC] = computeMeanTheta_resultantLength(theta,grpAxBxC);
    LRTS_AxBxC = kHat * sum(r_AxBxC .* (1 - abs(cos(theta_AxBxC - theta_grand))));
    R2_AxBxC = 1-exp(-LRTS_AxBxC / n);
end;

%--- test of significance of main and interaction effects-------
% LRTS (i.e. difference between log likelihood of 2 models) is approximately chi^2 distributed (Wilkin's theorem, see http://en.wikipedia.org/wiki/Likelihood-ratio_test)
df = 2 - 1; % in 2x2 design, we estimate 2 mean for alternative model and 1 grand mean for 0 model
p_A = 1-cdf('chi2',LRTS_A,df);
if nFact == 2 | nFact == 3;
    p_B = 1-cdf('chi2',LRTS_B,df);
    p_AxB = 1-cdf('chi2',LRTS_AxB,df);
end;
if nFact == 3;
    p_C = 1-cdf('chi2',LRTS_C,df);
    p_AxC = 1-cdf('chi2',LRTS_AxC,df);
    p_BxC = 1-cdf('chi2',LRTS_BxC,df);
    p_AxBxC = 1-cdf('chi2',LRTS_AxBxC,df);
end;


%--compile results table----
if nFact == 1;
    statTab{1,1} = grpNames{1};   
    statTab{1,2} = LRTS_A;    
    statTab{1,3} = p_A; 
    statTab{1,4} = R2_A;  
elseif nFact == 2;
    statTab{1,1} = grpNames{1};
    statTab{2,1} = grpNames{2};
    statTab{3,1} = [grpNames{1} ' X ' grpNames{2}] ;
    statTab{1,2} = LRTS_A;
    statTab{2,2} = LRTS_B;
    statTab{3,2} = LRTS_AxB;
    statTab{1,3} = p_A;
    statTab{2,3} = p_B;
    statTab{3,3} = p_AxB;
    statTab{1,4} = R2_A;
    statTab{2,4} = R2_B;
    statTab{3,4} = R2_AxB;      
elseif nFact == 3;
    statTab{1,1} = grpNames{1};
    statTab{2,1} = grpNames{2};
    statTab{3,1} = grpNames{3};
    statTab{4,1} = [grpNames{1} ' X ' grpNames{2}] ; % AxB
    statTab{5,1} = [grpNames{1} ' X ' grpNames{3}] ; % AxC
    statTab{6,1} = [grpNames{2} ' X ' grpNames{3}] ; % BxC
    statTab{7,1} = [grpNames{1} ' X ' grpNames{2} ' X ' grpNames{3}] ; % AxBxC
    statTab{1,2} = LRTS_A;
    statTab{2,2} = LRTS_B;
    statTab{3,2} = LRTS_C;
    statTab{4,2} = LRTS_AxB;
    statTab{5,2} = LRTS_AxC;
    statTab{6,2} = LRTS_BxC;
    statTab{7,2} = LRTS_AxBxC;
    statTab{1,3} = p_A;
    statTab{2,3} = p_B;
    statTab{3,3} = p_C;
    statTab{4,3} = p_AxB;
    statTab{5,3} = p_AxC;
    statTab{6,3} = p_BxC;
    statTab{7,3} = p_AxBxC;
    statTab{1,4} = R2_A;
    statTab{2,4} = R2_B;
    statTab{3,4} = R2_C;
    statTab{4,4} = R2_AxB;
    statTab{5,4} = R2_AxC;
    statTab{6,4} = R2_BxC;
    statTab{7,4} = R2_AxBxC;    
end;

function [meanTheta,resultantLength] = computeMeanTheta_resultantLength(theta,grp);
% computes mean theta and resultant length for different groups of data as
% indicated by grp. Important: grp must be 1,2,..k, k = number of groups
nGrps = length(unique(grp));
meanTheta = zeros(1,nGrps);
resultantLength = zeros(1,nGrps);
for iGrp = 1:nGrps;
    ind = grp == iGrp;
    s = sum(sin(theta(ind)));
    c = sum(cos(theta(ind)));
    if c > 0;
        meanTheta(iGrp) = atan(s/c);
    elseif (c == 0) && (s > 0);
        meanTheta(iGrp) = pi/2;
    elseif (c == 0) && (s < 0);
        meanTheta(iGrp) = -pi/2;
    else
        meanTheta(iGrp) = pi + atan(s/c);
    end;
    resultantLength(iGrp) = sqrt(c^2+s^2);
end;