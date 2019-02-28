function out = whos(obj)
% Method to keep the functionality of the Matlab.io.matFile object's whos method
% 
% USAGE:
%   out = whos(obj)
% INPUT:
%   obj (object): mvpares object
% OUTPUT: 
%   out (structure array): details of the fields of the dataset
%
out = whos(obj.data);
end