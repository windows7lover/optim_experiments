% load_function


paramFunction.x0 = zeros(nFeatures+1,1);

finfo = getFunction('logistic' , paramFunction);

% finfo.xstar = xstar;
% finfo.fstar = fstar;