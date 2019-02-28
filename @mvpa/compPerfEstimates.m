function varargout = compPerfEstimates(trueLab,predLab,mode)
% Method to compute various performance estimates across
% training and generalization samples.

% Checking if input data are matching
if (size(trueLab,1) ~= size(predLab,1))
    error('mvpa:compPerfEstimates:InputDimensionMismatch',...
        ['The number of examples in the ture and predicted ',...
        'label inputs must match!']);
end

[nExampl,nTrTimePoints,nGenTimePointsPerTrTimePoint] = size(predLab);

if strcmp(mode,'regression')
    
    [int,b,r2,r] = deal(NaN(nTrTimePoints,nGenTimePointsPerTrTimePoint));
    const = ones(nExampl,1);

    parfor i = 1:nTrTimePoints
        for j = 1:nGenTimePointsPerTrTimePoint
            actPredLabels = predLab(:,i,j);
            [est,~,~,~,stats] = regress(actPredLabels,[const,trueLab]);
            int(i,j) = est(1);
            b(i,j) = est(2);
            r2(i,j) = stats(1);
            r(i,j) = corr(trueLab,actPredLabels);
        end
    end

    varargout{1} = int;
    varargout{2} = b;
    varargout{3} = r2;
    varargout{4} = r;
    
elseif strcmp(mode,'classification')
    
    acc = NaN(nTrTimePoints,nGenTimePointsPerTrTimePoint);

    parfor i = 1:nTrTimePoints
        for j = 1:nGenTimePointsPerTrTimePoint
            actPredLabels = predLab(:,i,j);
            trueLabLevels = unique(trueLab);
            tempAcc = NaN(size(trueLabLevels));
            for k = 1:numel(trueLabLevels)
                idx = trueLab == trueLabLevels(k);
                tempAcc(k) = sum(trueLab(idx) == actPredLabels(idx))/sum(idx);
            end
            acc(i,j) = mean(tempAcc);
        end
    end

    varargout{1} = acc;
    
else
    error('mvpa:compPerfEstimates:invalidInput',...
          'The value of ''mode'' is invalid. ');
end

end