function avg = avgFTdata(varargin)
% Accessory function to compute the average of FT datasets
avg = varargin{1};
if isfield(avg,'individual')
    parameter = 'individual';
else
    parameter = avg.cfg.parameter{:};
end
temp = cellfun(@(x) x.(parameter),varargin,'UniformOutput',false);
dimToAvg = ndims(temp{1}) + 1;
avg.(parameter) = nanmean(cat(dimToAvg,temp{:}),dimToAvg);
end