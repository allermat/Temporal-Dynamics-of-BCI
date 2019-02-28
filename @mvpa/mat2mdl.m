function mdl = mat2mdl(M)
% Converts the libsvm's svmtrain output model converted to a matrix back to sturcture form. 
% 
% DETAILS:
%   CAUTION! 
%       - Only tested with libsvm-3.20  
%       - Probability estimates are not saved.
% 
% USAGE: 
%   mdl = mat2mdl(M)
%
% INPUT:
%   M: the data of libsvm's svmtrain output structure converted to a matrix
%       Dimensions of the matrix: 
%           nRows: 2 + (nClasses-1) + nFeatures 
%           nCols: number of totalSV or 8 or 12 depending on SV type and 
%                  number of totalSV
% 
% OUTPUT: 
%   mdl: libsvm's svmtrain output structure. 

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parse input
p = inputParser;

addRequired(p,'mdl',@ismatrix);
parse(p,M);

M = p.Results.mdl;

% svm model structure specification from libsvm3.20 README
% struct svm_model
% 	{
% 		struct svm_parameter param;	/* parameter */
% 		int nr_class;		/* number of classes, = 2 in regression/one class svm */
% 		int l;			/* total #SV */
% 		struct svm_node **SV;		/* SVs (SV[l]) */
% 		double **sv_coef;	/* coefficients for SVs in decision functions (sv_coef[k-1][l]) */
% 		double *rho;		/* constants in decision functions (rho[k*(k-1)/2]) */
% 		double *probA;		/* pairwise probability information */
% 		double *probB;
% 		int *sv_indices;    /* sv_indices[0,...,nSV-1] are values in [1,...,num_traning_data] to indicate SVs in the training set */

% 		/* for classification only */

% 		int *label;		/* label of each class (label[k]) */
% 		int *nSV;		/* number of SVs for each class (nSV[k]) */
% 					/* nSV[0] + nSV[1] + ... + nSV[k-1] = l */
% 		/* XXX */
% 		int free_sv;		/* 1 if svm_model is created by svm_load_model*/
% 					/* 0 if svm_model is created by svm_train */
% 	};

% Generating model structure
mdl = struct();
mdl.Parameters = M(1,1:5)';
mdl.nr_class = M(1,6);
mdl.totalSV = M(1,7);
lengthRho = mdl.nr_class*(mdl.nr_class-1)/2;
mdl.rho = M(1,7+(1:lengthRho))';
if M(1,1) > 1
    mdl.Label = [];
else
    mdl.Label = M(1,7+lengthRho+(1:mdl.nr_class))';
end
mdl.sv_indices = M(2,1:mdl.totalSV)';
mdl.ProbA = [];
mdl.ProbB = [];
if M(1,1) > 1
    mdl.nSV = [];
else
    mdl.nSV = M(1,7+lengthRho+mdl.nr_class+(1:mdl.nr_class))';
end
mdl.sv_coef = M(2+(1:mdl.nr_class-1),1:mdl.totalSV)';
mdl.SVs = sparse(M(2+mdl.nr_class:end,1:mdl.totalSV)');


end

