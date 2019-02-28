function out = who(obj)
% Method to keep the functionality of the Matlab.io.matFile object's who method
% 
% USAGE:
%   out = who(obj)
% INPUT:
%   obj (object): mvpares object
% OUTPUT: 
%   out (cell array): list of fields of the dataset
%
out = who(obj.data);
end