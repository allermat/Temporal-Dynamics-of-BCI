function c = chooseClosest(vals,x)
% Method for choosing array entry closest to number x
[~,idx] = min(abs(vals-x));
c = vals(idx);
end
