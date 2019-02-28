function [ciLow,ciHigh] = confmean(A,dim,alpha)
% Method for computing confidence interval of mean
% Based on this https://uk.mathworks.com/matlabcentral/answers/159417-how-to-calculate-the-confidence-interval
% Standard error of the mean
SEM = std(A,0,dim)/sqrt(size(A,dim));
% T-Score
ts = tinv(1-alpha,size(A,dim)-1);
% Confidence Intervals 
ciLow = mean(A,dim) - ts*SEM;
ciHigh = mean(A,dim) + ts*SEM;
end

