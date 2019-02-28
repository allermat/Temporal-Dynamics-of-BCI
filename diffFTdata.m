function diffx_y = diffFTdata(x,y)
% Accessory function to compute the difference of two FieldTrip datasets
if isfield(x,'individual')
    parameter = 'individual';
else
    parameter = x.cfg.parameter{:};
end
diffx_y = x;
diffx_y.(parameter) = x.(parameter) - y.(parameter);
end