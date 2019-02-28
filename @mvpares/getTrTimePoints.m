function trTimePoints = getTrTimePoints(obj)
% Method for getting the training time points
%
% Copyright(C) 2016, Mate Aller
if any(cellfun(@numel,obj.info.tr_timePoints) > 1)
    trTimePoints = cellfun(@mean,obj.info.tr_timePoints)';
else
    trTimePoints = cell2mat(obj.info.tr_timePoints)';
end
end