function setWritable(obj,isWritable)
% Method for setting the 'writable' property
% 
% USAGE:
%   setWritable(obj,isWritable)
% INPUT:
%   obj (object): mvpares object
%   isWritable (logical): whether to set the object writable
%
% Copyright(C) 2016, Mate Aller

% parsing input
p = inputParser;
addRequired(p,'obj');
addRequired(p,'isWritable',@(x)validateattributes(x,...
    {'logical'},{'scalar','nonempty'}));
parse(p,obj,isWritable);
obj = p.Results.obj;
isWritable = p.Results.isWritable;

obj.data.Properties.Writable = isWritable;
obj.writable = isWritable;

end