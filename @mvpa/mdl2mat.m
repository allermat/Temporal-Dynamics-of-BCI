function M = mdl2mat(mdl)
% Converts the libsvm's svmtrain output model to a matrix. 
% 
% DETAILS:
%   CAUTION! 
%       - Only tested with libsvm-3.20  
%       - Probability estimates are not saved. 
% 
% USAGE: 
%   M = mdl2mat(mdl)
%
% INPUT:
%   mdl: libsvm's svmtrain output structure. 
% 
% OUTPUT: 
%   M: the data of mdl converted to a matrix
%       Dimensions of the matrix: 
%           nRows: 2 + (nClasses-1) + nFeatures 
%           nCols: number of totalSV or 8 or 12 depending on SV type and 
%                  number of totalSV

% Copyright(C) 2016, Mate Aller
% allermat@gmail.com

% Parse input
p = inputParser;

fieldsReqMdl = {'Parameters','nr_class','totalSV','rho','Label','sv_indices',...
    'ProbA','ProbB','nSV','sv_coef','SVs'};
checkMdl = @(x) all(isfield(x,fieldsReqMdl));

addRequired(p,'mdl',checkMdl);
parse(p,mdl);

mdl = p.Results.mdl;

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

% If mdl.Parameters(1) > 1 mdl.Label and mdl.nSV is going to be empty
temp = [mdl.Parameters',mdl.nr_class,mdl.totalSV,mdl.rho',mdl.Label',mdl.nSV'];
% Matrix for storing mdl.Parameters, mdl.nr_class, mdl.totalSV, mdl.rho, 
% mdl.Label, mdl.nSV, mdl.sv_indices and mdl.sv_coef
% The number of columns in mdl.sv_coef is nr_class-1
M1 = NaN(2+mdl.nr_class-1,max([mdl.totalSV,size(temp,2)]));
M1(1,1:size(temp,2)) = temp;
M1(2,1:size(mdl.sv_indices,1)) = mdl.sv_indices';
M1(2+(1:mdl.nr_class-1),1:size(mdl.sv_coef,1)) = mdl.sv_coef';
% Matrix for storing mdl.SVs
M2 = NaN(size(mdl.SVs,2),max([mdl.totalSV,size(temp,2)]));
M2(:,1:mdl.totalSV) = full(mdl.SVs');
% Concatenating the two matrices
M = [M1;M2];

end

