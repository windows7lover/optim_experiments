function [ X,y,nFeatures,nPoints,w,fstar,reg ] = load_data( dataset, reg )
% dataset - reg : 
%   Madelon - 1e-3/1e2/1e7
%   Sido0   - 1e2
%   Sonar   - 1e-1/1e-6

dataname = [lower(dataset), '_data.mat'];
optvalname = [lower(dataset), '_opt_val_', num2str(reg), '.mat'];

load(dataname)
load(optvalname)



end
