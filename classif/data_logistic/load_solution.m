function [w,fstar] = load_solution( dataset, reg )
% dataset - reg : 
%   Madelon - 1e-3/1e2/1e7
%   Sido0   - 1e2
%   Sonar   - 1e-1/1e-6

optvalname = [lower(dataset), '_opt_val_', num2str(reg), '.mat'];
load(optvalname)



end