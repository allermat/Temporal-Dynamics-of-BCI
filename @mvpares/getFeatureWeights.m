function w = getFeatureWeights(obj)
% Method for accessing the trained weights for features
% 
% USAGE:
%   out = getFeatureWeights(obj)
% INPUT:
%   obj (object): mvpares object
% OUTPUT:
%   w (double): 3D array containing the feature weights for every
%       fold and timepoint. Size: nFeatures x nCVfolds x nTimePoints
%
% Copyright(C) 2018, Mate Aller

% Assign empty array and return if the result dataset is generalized
if strcmp(obj.state,'trained_and_generalized')
    % This extra check is for backwards compatibility. In older
    % versions tr_models were saved in the generalized files as
    % well. 
    if ~ismember('tr_models',obj.who)
        w = [];
        warning('mvpares:getFeatureWeights:datasetStateMismatch',...
                ['The dataset''s state is ''traned_and_generalized'',', ...
                 'so it does not contain the trained feature weights anymore.']);
        return;
    end
end

models = obj.data.tr_models;

w = NaN(size(models,1)-3,size(models,3),size(models,4));

for i = 1:size(models,3)
    
    for j = 1:size(models,4)
        mdl = mvpa.mat2mdl(models(:,:,i,j));
        w(:,i,j) = sum(full(mdl.SVs).*repmat(mdl.sv_coef,1,size(mdl.SVs,2)),1);
    end
    
end

end