function estimates = fitparams_pseu(genPredLabels,genExamples,smallDiscr,relV,ciAlpha)
% Method for fitting design parameters to predicted labels

[nExamples,nTrTimePoints,nGenTimePointsPerTrTimePoints] = size(genPredLabels);

genExamples.largeDiscr = abs(genExamples.locA-genExamples.locV) > smallDiscr;

r_d = genExamples.relV == relV(2) & ~genExamples.largeDiscr;
r_D = genExamples.relV == relV(2) & genExamples.largeDiscr;
R_d = genExamples.relV == relV(1) & ~genExamples.largeDiscr;
R_D = genExamples.relV == relV(1) & genExamples.largeDiscr;

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

designs = {designVR,designDisp};
designVarNames = {designVRvarNames,designDispVarNames};
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

