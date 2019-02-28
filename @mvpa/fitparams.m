function estimates = fitparams(genPredLabels,genExamples,factors,cParam)
% Method for fitting design parameters to predicted labels

% Processing input
smallDisp = cParam.smallDisp;
relV = cParam.relVlevels;
ciAlpha = cParam.ciAlpha;

[nExamples,nTrTimePoints,nGenTimePointsPerTrTimePoints] = size(genPredLabels);

genExamples.largeDisp = abs(genExamples.locA-genExamples.locV) > smallDisp;

% Logical vectors representing the levels of the factors
if all(ismember({'Disp','VR','Task'},factors))
    r_d_a = genExamples.relV == relV(2) & ~genExamples.largeDisp & genExamples.task == 1;
    r_d_v = genExamples.relV == relV(2) & ~genExamples.largeDisp & genExamples.task == 2;
    r_D_a = genExamples.relV == relV(2) & genExamples.largeDisp & genExamples.task == 1;
    r_D_v = genExamples.relV == relV(2) & genExamples.largeDisp & genExamples.task == 2;
    R_d_a = genExamples.relV == relV(1) & ~genExamples.largeDisp & genExamples.task == 1;
    R_d_v = genExamples.relV == relV(1) & ~genExamples.largeDisp & genExamples.task == 2;
    R_D_a = genExamples.relV == relV(1) & genExamples.largeDisp & genExamples.task == 1;
    R_D_v = genExamples.relV == relV(1) & genExamples.largeDisp & genExamples.task == 2;
    designMode = 1;
elseif all(ismember({'Disp','VR'},factors))
    r_d = genExamples.relV == relV(2) & ~genExamples.largeDisp;
    r_D = genExamples.relV == relV(2) & genExamples.largeDisp;
    R_d = genExamples.relV == relV(1) & ~genExamples.largeDisp;
    R_D = genExamples.relV == relV(1) & genExamples.largeDisp;
    designMode = 2;
elseif all(ismember({'Disp','Task'},factors))
    d_a = ~genExamples.largeDisp & genExamples.task == 1;
    d_v = ~genExamples.largeDisp & genExamples.task == 2;
    D_a = genExamples.largeDisp & genExamples.task == 1;
    D_v = genExamples.largeDisp & genExamples.task == 2;
    designMode = 3;
elseif ismember('Disp',factors)
    d = ~genExamples.largeDisp;
    D = genExamples.largeDisp;
    designMode = 4;
else
    error('mvpa:fitparams:invalidInput',...
        'The value of ''factors'' is invalid');
end

switch designMode
  case 1
    % Design for the 3-way interactions
    design3way = [ones(nExamples,1),zeros(nExamples,16)];
    design3way(r_d_a,2) = genExamples.locA(r_d_a); % A_r_d_a
    design3way(r_d_v,3) = genExamples.locA(r_d_v); % A_r_d_v
    design3way(r_D_a,4) = genExamples.locA(r_D_a); % A_r_D_a
    design3way(r_D_v,5) = genExamples.locA(r_D_v); % A_r_D_v
    design3way(R_d_a,6) = genExamples.locA(R_d_a); % A_R_d_a
    design3way(R_d_v,7) = genExamples.locA(R_d_v); % A_R_d_v
    design3way(R_D_a,8) = genExamples.locA(R_D_a); % A_R_D_a
    design3way(R_D_v,9) = genExamples.locA(R_D_v); % A_R_D_v
    design3way(r_d_a,10) = genExamples.locV(r_d_a); % V_r_d_a
    design3way(r_d_v,11) = genExamples.locV(r_d_v); % V_r_d_v
    design3way(r_D_a,12) = genExamples.locV(r_D_a); % V_r_D_a
    design3way(r_D_v,13) = genExamples.locV(r_D_v); % V_r_D_v
    design3way(R_d_a,14) = genExamples.locV(R_d_a); % V_R_d_a
    design3way(R_d_v,15) = genExamples.locV(R_d_v); % V_R_d_v
    design3way(R_D_a,16) = genExamples.locV(R_D_a); % V_R_D_a
    design3way(R_D_v,17) = genExamples.locV(R_D_v); % V_R_D_v
    design3wayVarNames = {...
        'A_r_d_a','A_r_d_v','A_r_D_a','A_r_D_v',...
        'A_R_d_a','A_R_d_v','A_R_D_a','A_R_D_v',...
        'V_r_d_a','V_r_d_v','V_r_D_a','V_r_D_v',...
        'V_R_d_a','V_R_d_v','V_R_D_a','V_R_D_v'};

    % Designs for the 2-way interactions
    % Design VR x Disparity
    designVRxDisp = [ones(nExamples,1),zeros(nExamples,8)];
    designVRxDisp(r_d_a | r_d_v,2) = genExamples.locA(r_d_a | r_d_v); % A_r_d
    designVRxDisp(r_D_a | r_D_v,3) = genExamples.locA(r_D_a | r_D_v); % A_r_D
    designVRxDisp(R_d_a | R_d_v,4) = genExamples.locA(R_d_a | R_d_v); % A_R_d
    designVRxDisp(R_D_a | R_D_v,5) = genExamples.locA(R_D_a | R_D_v); % A_R_D
    designVRxDisp(r_d_a | r_d_v,6) = genExamples.locV(r_d_a | r_d_v); % V_r_d
    designVRxDisp(r_D_a | r_D_v,7) = genExamples.locV(r_D_a | r_D_v); % V_r_D
    designVRxDisp(R_d_a | R_d_v,8) = genExamples.locV(R_d_a | R_d_v); % V_R_d
    designVRxDisp(R_D_a | R_D_v,9) = genExamples.locV(R_D_a | R_D_v); % V_R_D
    designVRxDispVarNames = {...
        'A_r_d','A_r_D','A_R_d','A_R_D',...
        'V_r_d','V_r_D','V_R_d','V_R_D'};

    % Design VR x Report
    designVRxRep = [ones(nExamples,1),zeros(nExamples,8)];
    designVRxRep(r_d_a | r_D_a,2) = genExamples.locA(r_d_a | r_D_a); % A_r_a
    designVRxRep(r_d_v | r_D_v,3) = genExamples.locA(r_d_v | r_D_v); % A_r_v
    designVRxRep(R_d_a | R_D_a,4) = genExamples.locA(R_d_a | R_D_a); % A_R_a
    designVRxRep(R_d_v | R_D_v,5) = genExamples.locA(R_d_v | R_D_v); % A_R_v
    designVRxRep(r_d_a | r_D_a,6) = genExamples.locV(r_d_a | r_D_a); % V_r_a
    designVRxRep(r_d_v | r_D_v,7) = genExamples.locV(r_d_v | r_D_v); % V_r_v
    designVRxRep(R_d_a | R_D_a,8) = genExamples.locV(R_d_a | R_D_a); % V_R_a
    designVRxRep(R_d_v | R_D_v,9) = genExamples.locV(R_d_v | R_D_v); % V_R_v
    designVRxRepVarNames = {...
        'A_r_a','A_r_v','A_R_a','A_R_v',...
        'V_r_a','V_r_v','V_R_a','V_R_v'};

    % Design Disparity x Report
    designDispXrep = [ones(nExamples,1),zeros(nExamples,8)];
    designDispXrep(r_d_a | R_d_a,2) = genExamples.locA(r_d_a | R_d_a); % A_d_a
    designDispXrep(r_d_v | R_d_v,3) = genExamples.locA(r_d_v | R_d_v); % A_d_v
    designDispXrep(r_D_a | R_D_a,4) = genExamples.locA(r_D_a | R_D_a); % A_D_a
    designDispXrep(r_D_v | R_D_v,5) = genExamples.locA(r_D_v | R_D_v); % A_D_v
    designDispXrep(r_d_a | R_d_a,6) = genExamples.locV(r_d_a | R_d_a); % V_d_a
    designDispXrep(r_d_v | R_d_v,7) = genExamples.locV(r_d_v | R_d_v); % V_d_v
    designDispXrep(r_D_a | R_D_a,8) = genExamples.locV(r_D_a | R_D_a); % V_D_a
    designDispXrep(r_D_v | R_D_v,9) = genExamples.locV(r_D_v | R_D_v); % V_D_v
    designDispXrepVarNames = {...
        'A_d_a','A_d_v','A_D_a','A_D_v',...
        'V_d_a','V_d_v','V_D_a','V_D_v'};

    % Designs for the main effects
    % Design VR
    designVR = [ones(nExamples,1),zeros(nExamples,4)];
    designVR(r_d_a | r_d_v | r_D_a | r_D_v,2) = ...
        genExamples.locA(r_d_a | r_d_v | r_D_a | r_D_v); % A_r
    designVR(R_d_a | R_d_v | R_D_a | R_D_v,3) = ...
        genExamples.locA(R_d_a | R_d_v | R_D_a | R_D_v); % A_R
    designVR(r_d_a | r_d_v | r_D_a | r_D_v,4) = ...
        genExamples.locV(r_d_a | r_d_v | r_D_a | r_D_v); % V_r
    designVR(R_d_a | R_d_v | R_D_a | R_D_v,5) = ...
        genExamples.locV(R_d_a | R_d_v | R_D_a | R_D_v); % V_R
    designVRvarNames = {'A_r','A_R','V_r','V_R'};

    % Design Disparity
    designDisp = [ones(nExamples,1),zeros(nExamples,4)];
    designDisp(r_d_a | r_d_v | R_d_a | R_d_v,2) = ...
        genExamples.locA(r_d_a | r_d_v | R_d_a | R_d_v); % A_d
    designDisp(r_D_a | r_D_v | R_D_a | R_D_v,3) = ...
        genExamples.locA(r_D_a | r_D_v | R_D_a | R_D_v); % A_D
    designDisp(r_d_a | r_d_v | R_d_a | R_d_v,4) = ...
        genExamples.locV(r_d_a | r_d_v | R_d_a | R_d_v); % V_d
    designDisp(r_D_a | r_D_v | R_D_a | R_D_v,5) = ...
        genExamples.locV(r_D_a | r_D_v | R_D_a | R_D_v); % V_D
    designDispVarNames = {'A_d','A_D','V_d','V_D'};

    % Design Report
    designRep = [ones(nExamples,1),zeros(nExamples,4)];
    designRep(r_d_a | r_D_a | R_d_a | R_D_a,2) = ...
        genExamples.locA(r_d_a | r_D_a | R_d_a | R_D_a); % A_a
    designRep(r_d_v | r_D_v | R_d_v | R_D_v,3) = ...
        genExamples.locA(r_d_v | r_D_v | R_d_v | R_D_v); % A_v
    designRep(r_d_a | r_D_a | R_d_a | R_D_a,4) = ...
        genExamples.locV(r_d_a | r_D_a | R_d_a | R_D_a); % V_a
    designRep(r_d_v | r_D_v | R_d_v | R_D_v,5) = ...
        genExamples.locV(r_d_v | r_D_v | R_d_v | R_D_v); % V_v
    designRepVarNames = {'A_a','A_v','V_a','V_v'};

    designs = {design3way,designVRxDisp,designVRxRep,designDispXrep,designVR,...
               designDisp,designRep};
    designVarNames = {design3wayVarNames,designVRxDispVarNames,designVRxRepVarNames,...
                      designDispXrepVarNames,designVRvarNames,designDispVarNames,designRepVarNames};    
  case 2
    % Designs for the main effects
    % Design VR
    designVR = [ones(nExamples,1),zeros(nExamples,4)];
    designVR(r_d | r_D,2) = genExamples.locA(r_d | r_D); % A_r
    designVR(R_d | R_D,3) = genExamples.locA(R_d | R_D); % A_R
    designVR(r_d | r_D,4) = genExamples.locV(r_d | r_D); % V_r
    designVR(R_d | R_D,5) = genExamples.locV(R_d | R_D); % V_R
    designVRvarNames = {'A_r','A_R','V_r','V_R'};

    % Design Disparity
    designDisp = [ones(nExamples,1),zeros(nExamples,4)];
    designDisp(r_d | R_d,2) = genExamples.locA(r_d | R_d); % A_d
    designDisp(r_D | R_D,3) = genExamples.locA(r_D | R_D); % A_D
    designDisp(r_d | R_d,4) = genExamples.locV(r_d | R_d); % V_d
    designDisp(r_D | R_D,5) = genExamples.locV(r_D | R_D); % V_D
    designDispVarNames = {'A_d','A_D','V_d','V_D'};
    
    % Designs for the 2-way interactions
    % Design VR x Disparity
    designVRxDisp = [ones(nExamples,1),zeros(nExamples,8)];
    designVRxDisp(r_d,2) = genExamples.locA(r_d); % A_r_d
    designVRxDisp(r_D,3) = genExamples.locA(r_D); % A_r_D
    designVRxDisp(R_d,4) = genExamples.locA(R_d); % A_R_d
    designVRxDisp(R_D,5) = genExamples.locA(R_D); % A_R_D
    designVRxDisp(r_d,6) = genExamples.locV(r_d); % V_r_d
    designVRxDisp(r_D,7) = genExamples.locV(r_D); % V_r_D
    designVRxDisp(R_d,8) = genExamples.locV(R_d); % V_R_d
    designVRxDisp(R_D,9) = genExamples.locV(R_D); % V_R_D
    designVRxDispVarNames = {...
        'A_r_d','A_r_D','A_R_d','A_R_D',...
        'V_r_d','V_r_D','V_R_d','V_R_D'};
    
    designs = {designVRxDisp,designVR,designDisp};
    designVarNames = {designVRxDispVarNames,designVRvarNames,designDispVarNames};
  case 3
    % Designs for the main effects
    % Design Disparity
    designDisp = [ones(nExamples,1),zeros(nExamples,4)];
    designDisp(d_a | d_v,2) = genExamples.locA(d_a | d_v); % A_d
    designDisp(D_a | D_v,3) = genExamples.locA(D_a | D_v); % A_D
    designDisp(d_a | d_v,4) = genExamples.locV(d_a | d_v); % V_d
    designDisp(D_a | D_v,5) = genExamples.locV(D_a | D_v); % V_D
    designDispVarNames = {'A_d','A_D','V_d','V_D'};

    % Design Report
    designRep = [ones(nExamples,1),zeros(nExamples,4)];
    designRep(d_a | D_a,2) = genExamples.locA(d_a | D_a); % A_a
    designRep(d_v | D_v,3) = genExamples.locA(d_v | D_v); % A_v
    designRep(d_a | D_a,4) = genExamples.locV(d_a | D_a); % V_a
    designRep(d_v | D_v,5) = genExamples.locV(d_v | D_v); % V_v
    designRepVarNames = {'A_a','A_v','V_a','V_v'};
    
    % Designs for the 2-way interactions
    % Design Disparity x Report
    designDispXrep = [ones(nExamples,1),zeros(nExamples,8)];
    designDispXrep(d_a,2) = genExamples.locA(d_a); % A_d_a
    designDispXrep(d_v,3) = genExamples.locA(d_v); % A_d_v
    designDispXrep(D_a,4) = genExamples.locA(D_a); % A_D_a
    designDispXrep(D_v,5) = genExamples.locA(D_v); % A_D_v
    designDispXrep(d_a,6) = genExamples.locV(d_a); % V_d_a
    designDispXrep(d_v,7) = genExamples.locV(d_v); % V_d_v
    designDispXrep(D_a,8) = genExamples.locV(D_a); % V_D_a
    designDispXrep(D_v,9) = genExamples.locV(D_v); % V_D_v
    designDispXrepVarNames = {...
        'A_d_a','A_d_v','A_D_a','A_D_v',...
        'V_d_a','V_d_v','V_D_a','V_D_v'};
    
    designs = {designDispXrep,designDisp,designRep};
    designVarNames = {designDispXrepVarNames,designDispVarNames,designRepVarNames};
  case 4
    % In this case there is just only one main effect for Disparity
    designDisp = [ones(nExamples,1),zeros(nExamples,4)];
    designDisp(d,2) = genExamples.locA(d); % A_d
    designDisp(D,3) = genExamples.locA(D); % A_D
    designDisp(d,4) = genExamples.locV(d); % V_d
    designDisp(D,5) = genExamples.locV(D); % V_D
    designDispVarNames = {'A_d','A_D','V_d','V_D'};
    designs = {designDisp};
    designVarNames = {designDispVarNames};
end

% Adding the 'overall' design which takes all A and V trials and
% computes one Wav
designOverall = [ones(nExamples,1),genExamples.locA, ...
                 genExamples.locV];
designOverallVarNames = {'A_oa','V_oa'};
designs = cat(2,designs,designOverall);
designVarNames = cat(2,designVarNames,{designOverallVarNames});

estimates = cell(size(designs));
fieldNames = cell(size(designs));
% Iterate through designs training time points and generailzation time
% points per training time points
for i = 1:numel(designs)
    actDesign = designs{i};
    actEstimates = NaN((size(actDesign,2)-1)*2,nTrTimePoints,nGenTimePointsPerTrTimePoints);
    parfor j = 1:nTrTimePoints
        for k = 1:nGenTimePointsPerTrTimePoints
            actPredLabels = genPredLabels(:,j,k);
            [b,bint] = regress(actPredLabels,actDesign,ciAlpha);
            % Ignoring intercept term
            b(1) = [];
            bint(1,:) = [];
            % Saving the parameter estimates and the margins of error
            actEstimates(:,j,k) = [b;bint(:,2)-b];
        end
    end
    if nGenTimePointsPerTrTimePoints > 1
        actEstimates = mat2cell(actEstimates,ones(size(actEstimates,1),1),nTrTimePoints,nGenTimePointsPerTrTimePoints);
        actEstimates = cellfun(@squeeze,actEstimates,'UniformOutput',false);
    else
        actEstimates = mat2cell(actEstimates,ones(size(actEstimates,1),1),nTrTimePoints);
        actEstimates = cellfun(@transpose,actEstimates,'UniformOutput',false);
    end
    estimates{i} = actEstimates';
    fieldNames{i} = [strcat(designVarNames{i},'_beta'),strcat(designVarNames{i},'_ci')];
end

fieldNames = cat(2,fieldNames{:});
estimates = cat(2,estimates{:});

estimates = cell2struct(estimates,fieldNames,2);

end

