function s = getSizeGenTime(obj)
% Method for getting the size of the generalization time points
% 
% USAGE:
%   s = getSizeGenTime(obj)
% INPUT:
%   obj (object): mvpares object
% OUTPUT:
%   s (vector): size of the generalization time points array
%
% Copyright(C) 2016, Mate Aller

s = size(obj.getGenTimePoints);
end