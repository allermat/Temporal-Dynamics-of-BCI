function [scaledFeats,varargout] = scalefeatures(feats,method,varargin)
% Scales the features according to the specified method. 
% 
% USAGE: 
%   [scaledFeats,scParamOut] = scalefeatures(feats,method)
%   scaledFeats = scalefeatures(feats,method,scParam)
% 
% INPUT: 
%   feats: NxM matrix of features. N = number of examples, M = number of
%       features.
%   method: 'Z-f','Z-e','Z-ef','Z-fe'. 
%           Z means Z-scoring
%           e: across examples (separately for features)
%           f: across features (separately for examples)
%   scParam: Nx2 matrix of parameters for Z scoring across examples. 
%            N = number of features
%            column 1: mean for Z-scoring
%            column 2: std for Z-scoring
% 
% OUTPUT: 
%   scaledFeats: the same set of features scaled accordig to the method. 
%   varargout: Nx2 matrix of parameters for Z scoring across examples
%              applied if they weren't specified in the input as scParam;
%              N = number of features
%              column 1: mean for Z-scoring
%              column 2: std for Z-scoring
%   

%% Parsing input
p = inputParser;

% Valid scaling methods
validScMethods = {'Z-f','Z-e','Z-ef','Z-fe'};

% Defining inputs.
addRequired(p,'feats',@(x)validateattributes(x,{'numeric'},{'nonempty'}));
addRequired(p,'method',@(x)any(validatestring(x,validScMethods)));
addOptional(p,'scParam',[],@(x)validateattributes(x,{'numeric'},{'size',[size(feats,2),2],'nonempty'}));

% Parsing inputs.
parse(p,feats,method,varargin{:});

% Assigning inputs to variables.
feats = p.Results.feats;
method = p.Results.method;
scParam = p.Results.scParam;

%% Main part of the function
scaledFeats = feats;
if isempty(scParam) && ~strcmp(method,'Z-f')
    scParamOut = NaN(size(feats,2),2);
end

if strcmp(method,'Z-f')
    
    % Normalizing across features (separately for examples)
    for i = 1:size(feats,1)
        scaledFeats(i,:) = (feats(i,:)-mean(feats(i,:)))/std(feats(i,:));
    end
    
elseif strcmp(method,'Z-e')
    
    % Normalizing across examples (separately for features)
    for i = 1:size(feats,2)
        if ~isempty(scParam)
            % Use the already specified scaling parameters for scaling
            scaledFeats(:,i) = (feats(:,i)-scParam(i,1))/scParam(i,2);
        else
            % Saving scaling parameters for later use and scale
            scParamOut(i,:) = [mean(feats(:,i)),std(feats(:,i))];
            scaledFeats(:,i) = (feats(:,i)-scParamOut(i,1))/scParamOut(i,2);
        end
    end
   
elseif strcmp(method,'Z-fe')
    
    % Normalizing across features (separately for examples)
    for i = 1:size(feats,1)
        scaledFeats(i,:) = (feats(i,:)-mean(feats(i,:)))/std(feats(i,:));
    end
    
    % Normalizing across examples (separately for features)
    for i = 1:size(scaledFeats,2)
        if ~isempty(scParam)
            scaledFeats(:,i) = (scaledFeats(:,i)-scParam(i,1))/scParam(i,2);
        else
            scParamOut(i,:) = [mean(scaledFeats(:,i)),std(scaledFeats(:,i))];
            scaledFeats(:,i) = (scaledFeats(:,i)-scParamOut(i,1))/scParamOut(i,2);
        end
    end
    
elseif strcmp(method,'Z-ef')
    
    % Normalizing across examples (separately for features)
    for i = 1:size(feats,2)
        if ~isempty(scParam)
            scaledFeats(:,i) = (feats(:,i)-scParam(i,1))/scParam(i,2);
        else
            scParamOut(i,:) = [mean(feats(:,i)),std(feats(:,i))];
            scaledFeats(:,i) = (feats(:,i)-scParamOut(i,1))/scParamOut(i,2);
        end
    end
    
    % Normalizing across features (separately for examples)
    for i = 1:size(scaledFeats,1)
        scaledFeats(i,:) = (scaledFeats(i,:)-mean(scaledFeats(i,:)))/std(scaledFeats(i,:));
    end
    
end

if isempty(scParam) && ~strcmp(method,'Z-f')
    % If there were no scaling parameters specified, output the ones used. 
    varargout{1} = scParamOut;
end


end