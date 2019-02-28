function fs = getFsample(obj)
% Method for accessing the training sampling frequency
% 
% USAGE:
%   fs = getFsample(obj)
% INPUT:
%   obj (object): mvpares object
% OUTPUT:
%   fs (scalar): estimated sampling fequency in Hz
%
% Copyright(C) 2016, Mate Aller

time = obj.getTrTimePoints;
% Tolerance for floating point comparison
tol = eps(0.5);
if any(abs(diff(diff(time))) > tol)
    warning('mvpares:getFsample:samplingNotUniform',...
        ['The sampling of the training time points is not uniform, the ',...
        'sampling frequency estimation might not be precise.']);
end
fs = 1/mean(diff(time));

end

