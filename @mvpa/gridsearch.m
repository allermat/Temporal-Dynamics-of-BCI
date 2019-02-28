function [bestparam,details] = gridsearch(TR_labs,TR_feats,svm_type,paramCoarse,nFold,nStepsFine)
% Performs grid search with two passes on the given parameters with nested cross-validation. 
% 
% USAGE: 
%   [bestparam,details] = gridsearch(TR_labs,TR_feats,svm_type,paramCoarse,nFold,nStepsFine)
% 
% INPUT: 
%   TR_labs: training labels (nExamples x 1)
%   TR_feats: training features (nExamples x nFeatures)
%   svm_type: type of SVM to be used ('er' for eSVR, 'nr' for nuSVR - these
%       two are supported at the moment with linear kernel). 
%   paramCoarse: cell array containing the parameter vectors for the course pass 
%       paramc{1,1} = param1 (nParam1 x 1)
%       paramc{1,2} = param2 (nParam2 x 1) (if applicable)
%   nFold: number of folds for the nested crossvalidation (scalar)
%   nStepsFine: number of steps for the fine pass of the grid search. 
%   
% OUTPUT: 
%   bestparam: vector with the best parameters overall (scalar or 2 x 1)
%   details: 2 x 4 cell array containing the details of the parameter
%       selection process separately for the two passes. 
%       First row - course pass, second row - fine pass.
%       {:,1} = CV MSE matrix nParam1 x nParam2 matrix
%       {:,2} = parameter 1 (nParam1 x 1)
%       {:,3} = parameter 2 (nParam2 x 1) (if applicable)
%       {:,4} = best paramters in the given pass (scalar or 2 x 1)
% 

if strcmp(svm_type,'cc')
    commArgs = sprintf('-s 0 -t 0 -v %d -q',nFold);
    param1flag = '-c';
    param1 = 2.^paramCoarse{1};
elseif strcmp(svm_type,'nc')
    commArgs = sprintf('-s 1 -t 0 -v %d -q',nFold);
    param1flag = '-n';
    param1 = paramCoarse{1};
elseif strcmp(svm_type,'er')
    commArgs = sprintf('-s 3 -t 0 -v %d -q',nFold);
    param1flag = '-c';
    param2flag = '-p';
    param1 = 2.^paramCoarse{1};
    param2 = 2.^paramCoarse{2};
elseif strcmp(svm_type,'nr')
    commArgs = sprintf('-s 4 -t 0 -v %d -q',nFold);
    param1flag = '-c';
    param2flag = '-n';
    param1 = 2.^paramCoarse{1};
    param2 = paramCoarse{2};
else
    error('The specified SVM type is not supported!');
end

% For support vector classification
if ismember(svm_type,{'cc','nc'})
    
    cvACCc = zeros(size(paramCoarse{1},2),1);
    
    nParam1 = size(param1,2);
        
    parfor i = 1:nParam1
        actparam1 = param1(i);
        args = sprintf('%s %s %g',commArgs,param1flag,actparam1);
        cvACCc(i) = svmtrain(TR_labs,TR_feats,args);
    end
    
    % Getting the linear index of the minimum MSE value.
    [~,idi] = max(cvACCc);
    % Determine best estimate in the coarse pass.
    bestparam = paramCoarse{1}(idi);
    % Saving the details
    details{1,1} = cvACCc;
    details{1,2} = paramCoarse{1}';
    details{1,3} = NaN;
    details{1,4} = bestparam;
    
    if idi == 1, paramFine{1} = paramCoarse{1}([idi idi+2]);
    elseif idi == length(paramCoarse{1}), paramFine{1} = paramCoarse{1}([idi-2 idi]);
    else paramFine{1} = paramCoarse{1}([idi-1 idi+1]);
    end
    
    paramFine{1} = paramFine{1}(1):((paramFine{1}(2)-paramFine{1}(1))/(nStepsFine-1)):paramFine{1}(2);
        
    cvACCf = zeros(size(paramFine{1},2),1);
    
    % Iterating through the parameter space, fine pass.
    if strcmp(svm_type,'cc')
        param1 = 2.^paramFine{1};
    elseif strcmp(svm_type,'nc')
        param1 = paramFine{1};
    end
    
    nParam1 = size(param1,2);
    
    parfor i = 1:nParam1
        actparam1 = param1(i);
        args = sprintf('%s %s %g',commArgs,param1flag,actparam1);
        cvACCf(i) = svmtrain(TR_labs,TR_feats,args);
    end
    
    % Getting the linear index of the minimum MSE value.
    [~,idi] = max(cvACCf);
    % Determine best estimate in the fine pass.
    bestparam = paramFine{1}(idi);
    % Saving the details
    details{2,1} = cvACCf;
    details{2,2} = paramFine{1}';
    details{2,3} = NaN;
    details{2,4} = bestparam;
    
% For support vector regression
else
    
    cvMSEc = zeros(size(paramCoarse{1},2),size(paramCoarse{2},2));
    
    nParam1 = size(param1,2);
    nParam2 = size(param2,2);
    
    for i = 1:nParam1
        actparam1 = param1(i);
        parfor j = 1:nParam2
            args = sprintf('%s %s %g %s %g',commArgs,param1flag,actparam1,param2flag,param2(j));
            cvMSEc(i,j) = svmtrain(TR_labs,TR_feats,args);
        end
    end
    
    % Getting the linear index of the minimum MSE value.
    [~,idl] = min(cvMSEc(:));
    % Computing the subscript indices according to the linear index.
    [idi,idj] = ind2sub(size(cvMSEc),idl);
    % Determine best estimate in the coarse pass.
    bestparam = [paramCoarse{1}(idi);paramCoarse{2}(idj)];
    % Saving the details
    details{1,1} = cvMSEc;
    details{1,2} = paramCoarse{1}';
    details{1,3} = paramCoarse{2}';
    details{1,4} = bestparam;
    
    if idi == 1, paramFine{1} = paramCoarse{1}([idi idi+2]);
    elseif idi == length(paramCoarse{1}), paramFine{1} = paramCoarse{1}([idi-2 idi]);
    else paramFine{1} = paramCoarse{1}([idi-1 idi+1]);
    end
    
    paramFine{1} = paramFine{1}(1):((paramFine{1}(2)-paramFine{1}(1))/(nStepsFine-1)):paramFine{1}(2);
    
    if idj == 1, paramFine{2} = paramCoarse{2}([idj idj+2]);
    elseif idj == length(paramCoarse{2}), paramFine{2} = paramCoarse{2}([idj-2 idj]);
    else paramFine{2} = paramCoarse{2}([idj-1 idj+1]);
    end
    
    paramFine{2} = paramFine{2}(1):((paramFine{2}(2)-paramFine{2}(1))/(nStepsFine-1)):paramFine{2}(2);
    
    cvMSEf = zeros(size(paramFine{1},2),size(paramFine{2},2));
    
    % Iterating through the combinations of the parameters in the parameter
    % space, fine pass.
    if strcmp(svm_type,'er')
        param1 = 2.^paramFine{1};
        param2 = 2.^paramFine{2};
    elseif strcmp(svm_type,'nr')
        param1 = 2.^paramFine{1};
        param2 = paramFine{2};
    end
    
    nParam1 = size(param1,2);
    nParam2 = size(param2,2);
    
    for i = 1:nParam1
        actparam1 = param1(i);
        parfor j = 1:nParam2
            args = sprintf('%s %s %g %s %g',commArgs,param1flag,actparam1,param2flag,param2(j));
            cvMSEf(i,j) = svmtrain(TR_labs,TR_feats,args);
        end
    end
    
    % Getting the linear index of the minimum MSE value.
    [~,idl] = min(cvMSEf(:));
    % Computing the subscript indices according to the linear index.
    [idi,idj] = ind2sub(size(cvMSEf),idl);
    % Determine best estimate in the fine pass.
    bestparam = [paramFine{1}(idi);paramFine{2}(idj)];
    % Saving the details
    details{2,1} = cvMSEf;
    details{2,2} = paramFine{1}';
    details{2,3} = paramFine{2}';
    details{2,4} = bestparam;
    
end

end