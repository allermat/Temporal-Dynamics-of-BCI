function posPix = posdeg2pix(posDeg,viewDist,res)
% Converts positions in visual angle to pixels given the viewing distance and resolution. 
% 
% USAGE: 
%   posPix = posdeg2pix(posDeg,viewDist,res) 
%
% INPUT:
%   posDeg: the position in visual angle,numeric. 
%   viewDist: The viewing distance in mm, i.e the distance of the 
%       observer's eyes from the monitor. (see DETAILS for assumptions),
%       scalar
%   res: screen resolution in pixels/mm, scalar
% 
% OUTPUT: 
%   posPix: the position in pixels. 
% 
% DETAILS:
%   It is assumed that the observer is point-like.

%% Checking Matlab version
if verLessThan('matlab', '8.1.0.604')
    error('In order to run the function, 2013a or newer version of Matlab is needed. ');
end

%% Parsing input
p = inputParser;

addRequired(p,'posDeg',@isnumeric);
addRequired(p,'viewDist',@isscalar);
addRequired(p,'res',@isscalar);

parse(p,posDeg,viewDist,res);

posDeg = p.Results.posDeg;
viewDist = p.Results.viewDist;
res = p.Results.res;

%%
posPix = (tand(posDeg)*viewDist)*res;

end

