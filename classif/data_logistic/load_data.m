function [ X,y,nFeatures,nPoints,Xtest,ytest,nPointsTest ] = load_data( dataset, precent_test, varargin )
% dataset - reg : 
%   Madelon - 1e-3/1e2/1e7
%   Sido0   - 1e2
%   Sonar   - 1e-1/1e-6

if(strcmpi(dataset,'random') )
    n = varargin{1};
    nFeatures = varargin{2};
    data = logistic_regression_data_generator(n, nFeatures+1);
    X = [data.x_train , data.x_test];
    y = [data.y_train' ; data.y_test'];
    nPoints = length(y);
else
    dataname = [lower(dataset), '_data_normalized.mat'];
    load(dataname)
end

perm_idx = randperm(nPoints);

nPointsTest = round(nPoints*precent_test);
nPoints = nPoints-nPointsTest;

Xtest = X(:,perm_idx(1:nPointsTest));
ytest = y(perm_idx(1:nPointsTest));


X = X(:,perm_idx(nPointsTest+1:end));
y = y(perm_idx(nPointsTest+1:end));



end