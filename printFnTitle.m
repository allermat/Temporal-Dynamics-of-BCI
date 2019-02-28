function printFnTitle(width,fnName,ds)
% Method for printing function output title on the console
nSpace = width-length(fnName)-length(ds);
fprintf('\n%s%s%s\n%s\n',fnName,repmat(' ',1,nSpace),ds,repmat('-',1,width));
end

